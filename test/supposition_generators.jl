using Supposition, Supposition.Data
using FLOYao
using Supposition.Data: Possibility

# TODO:
# - Rename all the matrix_generator, bitstring_generator etc to just matrices, bitstrings,...?

# ----------------------------------------
# Complex numbers, matrices and bitstrings
# ----------------------------------------
function ComplexFloats(::Type{T}=Float64; kwargs...) where {T<:AbstractFloat}
    return map(complex, Data.Floats{T}(; kwargs...), Data.Floats{T}(; kwargs...))
end

function matrix_generator(elements::Possibility, nrows::Int, ncols::Int)
    entries = Data.Vectors(elements; min_size=nrows*ncols, max_size=nrows*ncols)
    return map(f -> reshape(f, nrows, ncols), entries)
end

matrix_generator(el::Possibility, n::Int) = matrix_generator(el, n, n)

function sparse_matrix_generator(n::Int; kwargs...) 
    entries = Data.Vectors(Data.Just(0.0) | Data.Floats{Float64}(; kwargs...);
                           min_size=n^2, max_size=n^2)
    return map(flat -> FLOYao.SparseMatrixCOO(reshape(flat, n, n)), entries)
end

bitstring_generator(nq) = map(BitStr{nq}, Data.Integers(0, 2^nq - 1))
bitstring_generator(; max_nq=10) = Data.bind(bitstring_generator, Data.Integers(1, max_nq))

ordered_pair_generator(n::Int) = map(minmax, Data.Integers(1, n), Data.Integers(1, n))

# --------------------------------
# (number conserving) FLO circuits
# --------------------------------
"""
    flo_circuit_generator(nq::Int)

A suppositions Possibility for a FLO circuit on `nq` qubits
"""
function flo_circuit_generator(nq::Int)
    gate_gen = map(ordered_pair_generator(nq),
                   Data.Floats{Float64}(minimum=0.0, maximum=4π, nans=false),
                   Data.SampledFrom([:xx, :yy, :xy, :yx])
        ) do (i,j), θ, sym
        i == j && return rot(kron(nq, i => Z), θ)
        zs = [qb => Z for qb in i+1:j-1]
        generator = if sym == :xx
            kron(nq, i=>X, zs..., j=>X)
        elseif sym == :yy
            kron(nq, i=>Y, zs..., j=>Y)
        elseif sym == :xy
            kron(nq, i=>X, zs..., j=>Y)
        else
            kron(nq, i=>Y, zs..., j=>X) 
        end
        rot(generator, θ)
    end
    gates = Data.Vectors(gate_gen; min_size=0, max_size=2 * nq^2)
    map(gs -> foldl(push!, gs; init=chain(nq)), gates)
end

"""
    flo_circuit_generator(; max_nq::Int)

A suppositions Possibility for a FLO circuit on up to `max_nq` qubits
"""
function flo_circuit_generator(; max_nq::Int=10)
    Data.bind(flo_circuit_generator, Data.Integers(2, max_nq))
end

function nc_flo_circuit_generator(nq::Int)
    gate_gen = map(ordered_pair_generator(nq),
                   Data.Floats{Float64}(minimum=0.0, maximum=4π, nans=false),
                   Data.SampledFrom([:xx, :xy])
        ) do (i,j), θ, sym
        # TODO: I would like to use Rz here, but that gives a different global phase!
        i == j && return put(nq, i => shift(θ))
        zs = [qb => Z for qb in i+1:j-1]
        return if sym == :xx
            rot(kron(nq, i=>X, zs..., j=>X), θ) * rot(kron(nq, i=>Y, zs..., j=>Y), θ)
        else
            rot(kron(nq, i=>X, zs..., j=>Y), θ) * rot(kron(nq, i=>Y, zs..., j=>X), -θ)
        end
    end
    gates = Data.Vectors(gate_gen; min_size=0, max_size=2 * nq^2)
    map(gs -> foldl(push!, gs; init=chain(nq)), gates)
end

function nc_flo_circuit_generator(; max_nq::Int=10)
    Data.bind(nc_flo_circuit_generator, Data.Integers(2, max_nq))
end

# ---------------------------------
# (definite occupation) FLO states
# ---------------------------------
function flo_state_generator(T, nq::Integer)
    map(flo_circuit_generator(nq), bitstring_generator(nq)) do circuit, bits
        FLOYao.product_state(T, bits) |> circuit
    end
end
flo_state_generator(nq::Integer) = flo_state_generator(Float64, nq::Integer)

function flo_state_generator(T; max_nq=10)
    Data.bind(nq -> flo_state_generator(T, nq), Data.Integers(2, max_nq))
end
flo_state_generator(; max_nq=10) = flo_state_generator(Float64; max_nq=max_nq)

function nc_flo_state_generator(T, nq::Integer)
    map(nc_flo_circuit_generator(nq), bitstring_generator(nq)) do circuit, bits
        FLOYao.product_state(T, bits) |> circuit
    end
end
nc_flo_state_generator(nq::Integer) = nc_flo_state_generator(Float64, nq::Integer)

function nc_flo_state_generator(T; max_nq=10)
    Data.bind(nq -> nc_flo_state_generator(T, nq), Data.Integers(2, max_nq))
end
nc_flo_state_generator(; max_nq=10) = nc_flo_state_generator(Float64; max_nq=max_nq)

# ------------------------------------
# (number conserving) FLO hamiltonians
# ------------------------------------
function flo_hamiltonian_generator(T, nq)
    term_gen = map(ordered_pair_generator(nq),
                   Data.Floats{Float64}(minimum=-1e5, maximum=1e5, nans=false),
                   Data.SampledFrom([:xx, :yy, :xy, :yx])
        ) do (i,j), coeff, sym
        i == j && return coeff * kron(nq, i => Z)
        zs = [qb => Z for qb in i+1:j-1]
        term = if sym == :xx
            kron(nq, i=>X, zs..., j=>X)
        elseif sym == :yy
            kron(nq, i=>Y, zs..., j=>Y)
        elseif sym == :xy
            kron(nq, i=>X, zs..., j=>Y)
        else
            kron(nq, i=>Y, zs..., j=>X) 
        end
        coeff * term
    end
    terms = Data.Vectors(term_gen; min_size=1, max_size=2 * nq^2)
    map(ts -> Add(Vector{AbstractBlock{2}}(ts)), terms)
end
flo_hamiltonian_generator(nq::Integer) = flo_hamiltonian_generator(Float64, nq::Integer)

function flo_hamiltonian_generator(T; max_nq=10)
    Data.bind(nq -> flo_hamiltonian_generator(T, nq), Data.Integers(2, max_nq))
end
flo_hamiltonian_generator(; max_nq=10) = flo_hamiltonian_generator(Float64; max_nq=max_nq)


function nc_flo_hamiltonian_generator(T, nq)
    term_gen = map(ordered_pair_generator(nq),
                   Data.Floats{Float64}(minimum=-1e5, maximum=1e5, nans=false),
                   Data.SampledFrom([:xx, :xy])
        ) do (i,j), coeff, sym
        i == j && return coeff * kron(nq, i => Z)
        zs = [qb => Z for qb in i+1:j-1]
        return if sym == :xx
            coeff * (kron(nq, i=>X, zs..., j=>X) + kron(nq, i=>Y, zs..., j=>Y))
        else
            coeff * (kron(nq, i=>X, zs..., j=>Y) - kron(nq, i=>Y, zs..., j=>X))
        end
    end
    terms = Data.Vectors(term_gen; min_size=1, max_size=2 * nq^2)
    map(ts -> Add(Vector{AbstractBlock{2}}(ts)), terms)
end
nc_flo_hamiltonian_generator(nq::Integer) = nc_flo_hamiltonian_generator(Float64, nq::Integer)

function nc_flo_hamiltonian_generator(T; max_nq=10)
    Data.bind(nq -> nc_flo_hamiltonian_generator(T, nq), Data.Integers(2, max_nq))
end
nc_flo_hamiltonian_generator(; max_nq=10) = nc_flo_hamiltonian_generator(Float64; max_nq=max_nq)
