# ------------------------------
# The actual Register Definition
# ------------------------------
"""
    MajoranaReg{T} <: AbstractRegister{2}
    MajoranaReg(state::AbstractMatrix{T<:Real})

A register holding the "state" of a Majorana operators when propagating through
a FLO circuit as a `2n×2n` matrix.

#  Warning
The `MajoranaReg` constructor will not normalize the `state` matrix. It is 
recommended to use `FLOYao.zero_state` or `FLOYao.product_state` to produce
your initial state.
"""
mutable struct MajoranaReg{T<:Real} <: AbstractRegister{2}
    state::Matrix{T}
end


MajoranaReg(::Type{T}, n::Integer) where {T} = MajoranaReg(Matrix{T}(undef, 2n, 2n))
MajoranaReg(n::Integer) = MajoranaReg(Float64, n)
Yao.nqubits(reg::MajoranaReg) = size(reg.state, 1) ÷ 2
Yao.nqudits(reg::MajoranaReg) = size(reg.state, 1) ÷ 2
Yao.nactive(reg::MajoranaReg) = Yao.nqubits(reg)
Yao.nbatch(reg::MajoranaReg) = 1
Yao.nremain(reg::MajoranaReg) = 0
Base.copy(reg::MajoranaReg) = MajoranaReg(copy(reg.state))
Base.eltype(reg::MajoranaReg) = eltype(reg.state)
Yao.datatype(::MajoranaReg{T}) where {T} = T
Yao.state(reg::MajoranaReg) = reg.state

# The detailed version showing the contents of the register in e.g. 
# the jupyter cell output
function Base.show(io::IO, m::MIME"text/plain", reg::MajoranaReg)
    println(io, typeof(reg), " with $(Yao.nqubits(reg)) qubits:")
    show(io, m, reg.state)
end

# Less detailed version that is used e.g. in string interpolations
function Base.show(io::IO, reg::MajoranaReg)
    print(io, typeof(reg), "($(nqubits(reg)))")
end


"""
    majorana2arrayreg(reg::MajoranaReg)

Converts a `2n×2n` MajoranaReg `reg` into a `2^n` ArrayReg.

# Note
This implementation is not very clever and should mainly be used for debugging
purposes with small numbers of qubits. If ⟨Ω|U|Ω⟩ is close to zero, it is not
very accurate.
"""
function majorana2arrayreg(reg::MajoranaReg)
    nq = nqubits(reg)

    # praying here, that the zero_state has non-zero overlap with the 
    # new vacuum state.
    # This first bit gets areg into the new vacuum by piping it through
    # the projector U|Ω⟩⟨Ω|U^† 
    areg = uniform_state(Complex{eltype(reg)}, nq)
    for i in 1:nq
        γ_i1 = majoranaop(nq, reg.state[:,2i-1])
        γ_i2 = majoranaop(nq, reg.state[:,2i])
        circuit = 1im*chain(nq, γ_i2, γ_i1) + igate(nq)
        areg |> circuit
    end
    normalize!(areg)
    return areg
end

# ------------------------------
# Some convenience constructors
# TODO:
# - FLOYao.rand_state(bit_str)
# ------------------------------
"""
    zero_state([T=Float64,] n)

Create a Majorana register on `n` qubits in the vacuum state ``|Ω⟩`` with
storage type `T`.
"""
function zero_state(::Type{T}, n::Integer) where {T}
    return MajoranaReg(Matrix{T}(I(2n)))
end

zero_state(n) = zero_state(Float64, n)

function zero_state_like(reg::AbstractRegister)
    return zero_state(real(datatype(reg)), nqubits(reg))
end


"""
    product_state!(reg::MajoranaReg, bit_str::BitStr)

Put `reg` into the product state described by `bit_str`
"""
function product_state!(reg::MajoranaReg, bit_str::DitStr{2,N,IT}) where {N,IT}
    zero_state!(reg)
    for (qb, bit) in enumerate(bit_str)
        if Bool(bit)
            reg |> put(N, qb => X)
        end
    end
    return reg
end

"""
    product_state([T=Float64,] bit_str)

A Majorana register in the computational basis state specified by the 
`bit_str`.
"""
function product_state(::Type{T}, bit_str::DitStr{2,N,IT}) where {N,T,IT}
    reg = MajoranaReg(T, N)
    product_state!(reg, bit_str)
    return reg
end

product_state(bit_str) = product_state(Float64, bit_str)

"""
    zero_state!(reg::MajoranaReg)

Put `reg` into the computational zero state
"""
function zero_state!(reg::MajoranaReg{T}) where {T}
    fill!(reg.state, zero(T))
    reg.state[diagind(reg.state)] .= one(T)
    return reg
end

"""
    one_state!(reg::MajoranaReg)

Put `reg` into the all ones state
"""
function one_state!(reg::MajoranaReg{T}) where {T}
    fill!(reg.state, zero(T))
    for (i, j) in enumerate(diagind(reg.state))
        reg.state[j] = (i % 4 == 1 ||  i % 4 == 0) ? one(T) : -one(T)
    end
    return reg
end
