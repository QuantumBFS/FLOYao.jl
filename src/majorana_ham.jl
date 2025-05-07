#=
#  Authors:   Jan Lukas Bosse
#  Copyright: 2022 Phasecraft Ltd.
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
=#

# -----------------------------------
# Hamiltonians for expectation values
# -----------------------------------
"""
    MajoranaTerm{T}
    MajoranaTerm(coeff::Number, indices::Vector{Int})

A Majorana monomial represented by its coefficient and indices on
which it acts non-trivially.
"""
struct MajoranaTerm{T<:Number}
    coeff::T
    indices::Vector{Int}
end

Base.length(mt::MajoranaTerm) = length(mt.indices)

Base.:(*)(s::Number, term::MajoranaTerm) = MajoranaTerm(s * term.coeff, term.indices)


"""
    MajoranaSum{T}
    MajoranaSum(terms::Vector{MajoranaTerm{T}})
    MajoranaSum(terms::Dict{Vector{Int}, T})

A sum of `MajoranaTerm`s.

Use [`kron2majoranasum`](@ref) to construct it from Yao Blocks.
"""
struct MajoranaSum{T<:Number}
    terms::Dict{Vector{Int},T}

    function MajoranaSum(terms::Dict{Vector{Int},T}) where {T}
        new{T}(terms)
    end

    function MajoranaSum(terms::AbstractVector{MajoranaTerm{T}}) where {T}
        firstterm, otherterms = Iterators.peel(terms)
        out = new{T}(Dict(firstterm.indices => firstterm.coeff))
        for term in otherterms
            addterm!(out, term)
        end
        return out
    end

    global function unsafe_majoranasum(terms::AbstractVector{MajoranaTerm{T}}) where {T}
        return new{T}(Dict(term.indices => term.coeff for term in terms))
    end
end

Base.copy(ms::MajoranaSum) = MajoranaSum(copy(ms.terms))

function addterm!(ms::MajoranaSum{T}, term::MajoranaTerm{T}) where {T}
    term.coeff ≈ 0 && return ms

    if term.indices in keys(ms.terms)
        new_coeff = ms.terms[term.indices] + term.coeff
        if new_coeff ≈ 0
            delete!(ms.terms, term.indices)
        else
            ms.terms[term.indices] = new_coeff
        end
    else
        ms.terms[term.indices] = term.coeff
    end

    return ms
end

function Base.iterate(ms::MajoranaSum, state...)
    inds_coeff_state = Base.iterate(ms.terms, state...)
    isnothing(inds_coeff_state) && return nothing
    ((indices, coeff), state) = inds_coeff_state
    return MajoranaTerm(coeff, indices), state
end

Base.length(ms::MajoranaSum) = length(ms.terms)

Base.eltype(::MajoranaSum{T}) where {T} = MajoranaTerm{T}
Base.eltype(::Type{MajoranaSum{T}}) where {T} = MajoranaTerm{T}

function Base.show(io::IO, ms::MajoranaSum)
    buf = IOBuffer()
    nterms = length(ms)
    println(buf, "$nterms-element $(typeof(ms)):")
    for mt in ms
        println(buf, "  ", mt.indices, " : ", mt.coeff)
    end
    print(io, String(take!(buf))[1:end-1]) # remove last newline
end

Base.zero(::Type{MajoranaSum{T}}) where {T} = MajoranaSum(Dict{Vector{Int}, T}())
Base.zero(::MajoranaSum{T}) where {T} = MajoranaSum(Dict{Vector{Int}, T}())

function _add_majoranasums(summands::NTuple{N,MajoranaSum{T}}) where {N,T}
    firstsum, othersums = Iterators.peel(summands)
    out = copy(firstsum)
    for mt in Iterators.flatten(othersums)
        addterm!(out, mt)
    end
    return out
end
Base.:(+)(summands::Vararg{MajoranaSum}) = _add_majoranasums(promote(summands...))
Base.:(+)(ms1::MajoranaSum, ms2::MajoranaSum) = _add_majoranasums(promote(ms1, ms2))
Base.:(+)(terms::Vararg{MajoranaTerm}) = MajoranaSum([terms...])

Base.:(*)(n::Number, ms::MajoranaSum) = unsafe_majoranasum(n .* ms)

"""
    kron2majoranaterm(::Type{T}, k::KronBlock)

Convert a `PauliKronBlock` to a `MajoranaTerm{complex(T)}`
"""
function kron2majoranaterm(::Type{T}, k::PauliKronBlock) where {T}

    locs = first.(k.locs)
    perm = YaoBlocks.TupleTools.sortperm(locs)

    # Stores the Majorana indices on which `k` acts non-trivially
    indices = Int[]
    sign = one(complex(T))
    # Stores if there is an even number of MajoranaOps after the current qubit
    outside_string = true

    # early out for identity
    length(locs) == 0 && return MajoranaTerm(sign, indices)

    q0, other_qubits = Iterators.peel(Iterators.reverse(perm))
    op = k.blocks[q0]
    if op == X
        push!(indices, 2locs[q0] - 1)
        outside_string = false
    elseif op == Y
        push!(indices, 2locs[q0])
        sign *= -1
        outside_string = false
    elseif op == Z
        push!(indices, 2locs[q0])
        push!(indices, 2locs[q0] - 1)
        sign *= -1im
    elseif op == I2
    else
        throw(DomainError(op, "Can only handle tensor products of Pauli operators as Hamiltonians"))
    end

    for q in other_qubits
        op = k.blocks[q]
        if outside_string
            if op == X
                push!(indices, 2locs[q] - 1)
                outside_string = false
            elseif op == Y
                push!(indices, 2locs[q])
                sign *= -1
                outside_string = false
            elseif op == Z
                push!(indices, 2locs[q])
                push!(indices, 2locs[q] - 1)
                sign *= -1im
            elseif op == I2
            else
                throw(DomainError(op, "Can only handle tensor products of Pauli operators as Hamiltonians"))
            end
        else
            # Fill in majorana string
            append!(indices, 2locs[q0]-2:-1:2locs[q]+1)
            sign *= (1im)^(locs[q0] - locs[q] - 1)
            if op == X
                push!(indices, 2locs[q])
                sign *= 1im
                outside_string = true
            elseif op == Y
                push!(indices, 2locs[q] - 1)
                sign *= 1im
                outside_string = true
            elseif op == Z
                sign *= -1
            elseif op == I2
                push!(indices, 2locs[q])
                push!(indices, 2locs[q] - 1)
            else
                throw(DomainError(op, "Can only handle tensor products of Pauli operators as Hamiltonians"))
            end
        end
        q0 = q
    end

    if !outside_string # shouldn't happen for even parity observables
        append!(indices, 2locs[q0]-2:-1:1)
        sign *= (-1im)^(locs[q0] - 1)
    end
    return MajoranaTerm(sign, reverse(indices))
end

"""
    yaoham2majoranasum(::Type{T}=Float64, yaoham::AbstracBlock{2})

Convert a  YaoBlock into a `MajoranaSum{complex(T)}`
"""
function yaoblock2majoranasum(::Type{T}, yaoham::Add{2}) where {T<:Real}
    return sum(yaoham, init=zero(MajoranaSum{complex(T)})) do block
        yaoblock2majoranasum(T, block)
    end
end

function yaoblock2majoranasum(::Type{T}, yaoham::KronBlock{2}) where {T<:Real}
    return unsafe_majoranasum([kron2majoranaterm(T, yaoham),])
end

function yaoblock2majoranasum(::Type{T}, yaoham::PutBlock{2,1,ZGate}) where {T<:Real}
    qb = yaoham.locs[1]
    return unsafe_majoranasum([MajoranaTerm(-1im * one(T), [2qb-1, 2qb])])
end

function yaoblock2majoranasum(::Type{T}, yaoham::PutBlock{2,N,<:PauliKronBlock}) where {T<:Real, N}
    mt = kron2majoranaterm(T, yaoham.content)
    # TODO: Not sure this logic is completely correct
    mt.indices .+= 2(minimum(yaoham.locs) - 1)
    return unsafe_majoranasum([mt])
end

function yaoblock2majoranasum(::Type{T}, yaoham::Scale) where {T<:Real}
    return factor(yaoham) * yaoblock2majoranasum(T, yaoham.content)
end

yaoblock2majoranasum(yaoham) = yaoblock2majoranasum(Float64, yaoham)


# -------------------------------
# Hamiltonians for time evolution
# -------------------------------
"""
    kron2majoranaindices(k::KronBlock)

Get the two(!) indices of Majorana operators from a kronecker product of pauli operators
"""
function kron2majoranasquare(k::KronBlock{2, M, <:NTuple{M, PauliGate}}) where {M}
    locs = first.(k.locs)
    perm = YaoBlocks.TupleTools.sortperm(locs)
    areconsecutive(ntuple(i->locs[perm[i]], M)) || throw(NonFLOException("$k is not acting on consecutive qubits"))

    # these are the indices in the majorana hamiltonian corresponding to this
    # kronecker product. swap accumulates the overall sign
    firstindex, secondindex, swap = -1, -1, true 
    for q in perm
        op = k.blocks[q]
        if firstindex == -1
            if op == X
                firstindex = 2locs[q]
            elseif op == Y
                firstindex = 2locs[q] - 1
            elseif op == Z
                firstindex = 2locs[q]
                secondindex = 2locs[q]-1
            elseif op == I2
            else
                throw(NonFLOException("$(k.blocks) acting on $(k.locs) is not FLO"))
            end
        elseif secondindex == -1
            if op == X
                secondindex = 2locs[q] - 1
            elseif op == Y
                secondindex = 2locs[q]
                swap = !swap
            elseif op == Z
                swap = !swap
            elseif op == I2
            else
                throw(NonFLOException("$(k.blocks) acting on $(k.locs) is not FLO"))
            end
        elseif op != I2
            throw(NonFLOException("$(k.blocks) acting on $(k.locs) is not FLO"))
        end
    end
    firstindex != -1 && secondindex != -1 || throw(NonFLOException("$(k.blocks) acting on $(k.locs) is not FLO"))
    # swap if we picked up a minus sign in total
    return swap ? (secondindex, firstindex) : (firstindex, secondindex)
end

kron2majoranasquare(k::PutBlock{2,1,ZGate}) = (2k.locs[1]-1, 2k.locs[1])

# TODO: Check if it is faster to build a `MajoranaSum` first and then convert
# it to a sparse 2n×2n matrix

"""
    yaoham2majoranasquares(::Type{T}=Float64, yaoham::AbstracBlock{2})

Convert a hamiltonian written as a YaoBlock into the corresponding 
``2n×2n`` majorana hamiltonian.
"""
function yaoham2majoranasquares(::Type{T}, yaoham::Add{2}) where {T<:Real}
    ham = zeros(T, 2nqubits(yaoham), 2nqubits(yaoham))
    @inbounds for k in yaoham
        if k isa Scale
            fast_add!(ham, _rmul!(yaoham2majoranasquares(T, k.content), k.alpha))
            #ham += k.alpha * yaoham2majoranasquares(T, k.content)
        elseif k isa KronBlock
            i1, i2 = kron2majoranasquare(k)
            ham[i1,i2] += 2
            ham[i2,i1] -= 2
        elseif k isa Add
            fast_add!(ham, yaoham2majoranasquares(T, k))
        elseif k isa PutBlock{2,1,ZGate}
            i1, i2 = kron2majoranasquare(k)
            ham[i1,i2] += 2
            ham[i2,i1] -= 2
        elseif k isa PutBlock{2,<:Any,<:PauliKronBlock}
            areconsecutive(k.locs) || throw(NonFLOException("$(k.content) contains terms that are not the product of two Majoranas"))
            i1, i2 = 2 * (minimum(k.locs) - 1) .+ kron2majoranasquare(k.content)
            ham[i1,i2] += 2
            ham[i2,i1] -= 2
        else
            throw(NonFLOException("$k not recognized as a square of Majorana operators"))
        end
    end
    return ham
end

function yaoham2majoranasquares(::Type{T}, yaoham::KronBlock{2}) where {T<:Real}
    i1, i2 = kron2majoranasquare(yaoham)
    return SparseMatrixCOO([i1, i2], [i2, i1], T[2, -2], 2nqubits(yaoham), 2nqubits(yaoham))
end

function yaoham2majoranasquares(::Type{T}, yaoham::PutBlock{2,1,ZGate}) where {T<:Real}
    i1, i2 = 2yaoham.locs[1]-1, 2yaoham.locs[1]
    return SparseMatrixCOO([i1, i2], [i2, i1], T[2, -2], 2nqubits(yaoham), 2nqubits(yaoham))
end

function yaoham2majoranasquares(::Type{T}, yaoham::PutBlock{2,N,<:PauliKronBlock}) where {T<:Real, N}
    areconsecutive(yaoham.locs) || throw(NonFLOException("$(yaoham.content) contains terms that are not the product of two Majoranas"))
    i1, i2 = 2 * (minimum(yaoham.locs) - 1) .+ kron2majoranasquare(yaoham.content)
    return SparseMatrixCOO([i1, i2], [i2, i1], T[2, -2], 2nqubits(yaoham), 2nqubits(yaoham))
end

function yaoham2majoranasquares(::Type{T}, yaoham::Scale) where {T<:Real}
    return yaoham.alpha * yaoham2majoranasquares(T, yaoham.content)
end

yaoham2majoranasquares(yaoham) = yaoham2majoranasquares(Float64, yaoham)
