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

# ------------------------------
# The actual Register Definition
# ------------------------------
"""
    MajoranaReg{T} <: AbstractRegister{2}
    MajoranaReg(state::AbstractMatrix{T<:Real})

A register holding the "state" of a Majorana operators when propagating through
a FLO circuit as a `2n×2n` matrix.

# Warning
The `MajoranaReg` constructor will not initialize the `state` matrix. It is 
recommended to use `FLOYao.zero_state` or `FLOYao.product_state` to produce
your initial state.
"""
struct MajoranaReg{T<:Real} <: AbstractRegister{2}
    state::Matrix{T}
end


MajoranaReg(::Type{T}, n::Integer) where {T} = MajoranaReg(Matrix{T}(undef, 2n, 2n))
MajoranaReg(n::Integer) = MajoranaReg(Float64, n)
Yao.nqubits(reg::MajoranaReg) = size(reg.state, 1) ÷ 2
Yao.nqudits(reg::MajoranaReg) = size(reg.state, 1) ÷ 2
Yao.nactive(reg::MajoranaReg) = Yao.nqubits(reg)
Yao.nbatch(reg::MajoranaReg) = 1
Yao.nremain(reg::MajoranaReg) = 0
Yao.state(reg::MajoranaReg) = reg.state
Base.eltype(::MajoranaReg{T}) where {T} = T
Yao.datatype(::MajoranaReg{T}) where {T} = T
Base.copy(reg::MajoranaReg) = MajoranaReg(copy(reg.state))
Base.similar(reg::MajoranaReg) = MajoranaReg(similar(reg.state))

function Base.copyto!(dst::MajoranaReg, src::MajoranaReg)
    nqubits(dst) != nqubits(src) && throw(DimensionMismatch("nqubits(dst) = $(nqubits(dst)) != nqubits(src) = $(nqubits(src))"))
    copyto!(state(dst), state(src))
end


function Base.:(==)(lhs::MajoranaReg, rhs::MajoranaReg)
    return nqubits(lhs) == nqubits(rhs) && state(lhs) == state(rhs)
end

# Only a sensible operator for Tangent registers. Create a separate type for 
# these to ensure I can not multiply primal registers with a scalar?
function Base.:(*)(number::Real, reg::MajoranaReg)
    return MajoranaReg(number .* reg.state)
end

function Base.isapprox(lhs::MajoranaReg, rhs::MajoranaReg)
    return nqubits(lhs) == nqubits(rhs) && isapprox(state(lhs), state(rhs))
end

# The detailed version showing the contents of the register in e.g. 
# the jupyter cell output
function Base.show(io::IO, ::MIME"text/plain", reg::MajoranaReg)
    println(io, typeof(reg), " with $(Yao.nqubits(reg)) qubits:")
    Base.print_array(io, state(reg))
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
purposes with small numbers of qubits. It pipes a random state ``|ψ⟩`` 
through the projector ``U|Ω⟩⟨Ω|U^†`` which may give inaccurate results if 
``⟨ψ|U|ψ⟩`` is very small.
"""
function majorana2arrayreg(reg::MajoranaReg)
    nq = nqubits(reg)

    # praying here, that the rand_state has non-zero overlap with the 
    # new vacuum state.
    # This first bit gets areg into the new vacuum by piping it through
    # the projector U|Ω⟩⟨Ω|U^† 
    areg = Yao.rand_state(Complex{eltype(reg)}, nq)
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

"""
    zero_state_like(reg::AbstractRegister)

Create a Majorana register in the zero state with the same element type 
and number of qubits as `reg`.
"""
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
    product_state([T=Float64,] bit_str::DitStr{2})
    product_state([T=Float64,] bit_configs::AbstractVector)
    product_state([T=Float64,] nbits::Int, val::Int)

Create an `MajoranaReg` of a product state.

The state can be specified as a bit string, as an array of Integers or Booleans
or with `nbits` and `val`.
"""
function product_state(::Type{T}, bit_str::DitStr{2,N,IT}) where {N,T,IT}
    reg = MajoranaReg(T, N)
    product_state!(reg, bit_str)
    return reg
end

product_state(bit_str) = product_state(Float64, bit_str)

# vector input
function product_state(::Type{T}, bit_configs::AbstractVector) where {T}
    # have to do conversion to DitStr ourselves because of integer overflow
    bit_str = sum(pairs(bit_configs), init=zero(BigInt)) do (i, b)
        b * BigInt(2)^(i-1)
    end |> DitStr{2,length(bit_configs),BigInt}
    return product_state(T, bit_str)
end

# integer input
function product_state(::Type{T}, nbits::Int, val::Integer) where {T}
    return product_state(T, DitStr{2,nbits}(val))
end


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
    rand_state([Float64,] n)

Create a Haar random `MajoranaReg` on `n` qubits.
"""
function rand_state(::Type{T}, n) where {T}
    state = random_orthogonal_matrix(T, 2n)
    return MajoranaReg(state)
end

rand_state(n) = rand_state(Float64, n)

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

"""
    one_state([T=Float64,] n)

Create a Majorana register on `n` qubits in the all one state ``|1 ⋯ 1⟩`` with
storage type `T`.
"""
function one_state(::Type{T}, n::Integer) where {T}
    reg = MajoranaReg(T, n)
    one_state!(reg)
    return reg
end

one_state(n) = one_state(Float64, n)
