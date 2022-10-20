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

# -----------------
# General utilities
# -----------------
const σs = [[1 0; 0 1], [0 1; 1 0], [0 -1im; 1im 0], [1 0; 0 -1]]
⊗ = kron

"""Tests if all entries in `locs` are consecutive"""
areconsecutive(locs) = all(locs[i+1] == locs[i]+1 || locs[i+1] == locs[i]-1 for i in 1:length(locs)-1)

struct NonFLOException <: Exception
    msg::String
end

Base.showerror(io::IO, exc::NonFLOException) = print(io, "NonFLOException: ", exc.msg)

"Majorana operator in a system with `n` on site Majorana site `i`"
function majoranaop(n, i::Integer)
    if i % 2 == 1 # so this one is of the first type (also: fuck 1-based indexing)
        return (-1)^((i-1)÷2) * kron(n, (j => Z for j in 1:(i-1)÷2)..., ((i+1)÷2 => X))
    else
        return -(-1)^((i-1)÷2) * kron(n, (j => Z for j in 1:(i-1)÷2)..., ((i+1)÷2 => Y))
    end
end

"Majorana op on `n` Majorana sites for orbital `ψ`"
function majoranaop(n, ψ::AbstractVector)
    return sum([ψ[i] * majoranaop(n, i) for i in eachindex(ψ)])
end

# -------------------------------------------------------------------
# Utilities to change between different bases for hermitian operators
# -------------------------------------------------------------------
@doc raw"""
    qubit2paulibasis(A::AbstractMatrix)

Converts a ``2^n×2^n`` matrix `A` in the standard qubit basis into a 
``4^n`` vector representing the same operator in the Pauli basis.

The ordering is as follows: Let ``σ^0 = I, σ^1 = X, σ^2 = Y`` and ``σ^3 = Z``. 

# Todo
Creating the dense pauli tensor product via `kron` and inner product via `dot`
is much slower than neccessary, since the pauli tensor product matrix is 
very sparse. Much faster would be to compute the inner product directly.
"""
function qubit2paulibasis(A::AbstractMatrix{T}) where {T}
    n = Int(log2(size(A, 1)))
    out = Vector{complex(T)}(undef, 4^n)
    for i in 0x0:UInt32(4^n - 1)
        # creating the tensor product of pauli matrices associated to i
        σ = mapreduce(k -> σs[((i >> 2k) % 4) + 1], kron, 0x0:UInt8(n-1), init=[one(T)])
        out[i+1] = σ ⋅ A / 2^n
    end
    return out
end

"""
    paulibasis2qubitop(P::AbstractVector)

Converts an operator to in the pauli basis to a matrix in
the computational basis. Inverse to `qubit2paulibasis`.
"""
function paulibasis2qubitop(P::AbstractVector{T}) where {T}
    n = Int(log(4, length(P)))
    out = zeros(complex(T), 2^n, 2^n)
    for i in 0x0:UInt32(4^n-1)
        σ = mapreduce(k -> σs[((i >> 2k) % 4) + 1], kron, 0x0:UInt8(n-1), init=[one(T)])
        out .+= P[i+1] .* σ
    end
    return out
end

"""
    paulibasis2majoranasquares(P::AbstractVector, locs=1:log4(length(P)))

Convert an operator written in the Pauli basis as a ``4^n``-element vector
to the corresponding ``2n×2n`` matrix of coefficients of products of two
Majorana operators.

Throws a `NonFLOException` if `P` contains terms that are not corresponding to 
the product of two Majorana operators.
"""
function paulibasis2majoranasquares(P::AbstractVector, locs=1:Int(log(4, length(P))))
    n = Int(log(4, length(P)))
    abs(P[1]) >= √(eps(abs(P[1]))) && @debug "Ignoring parts of the Hamiltonian that are proportional to the identity"
    
    areconsecutive(locs) || throw(NonFLOException("P contains terms that are not the product of two Majoranas"))
    
    out = zeros(real(eltype(P)), 2n, 2n)
    for i in 0x1:UInt32(4^n-1)
        coeff = P[i+1]
        if abs(coeff) >= √eps(eltype(out)) # really only look at non-zero coeffs
            firstmajorana = -1 # these will index into the 2n×2n matrix holding the majorana coeffs
            secondmajorana = -1
            for k in n-1:-1:0 
                opi = i % 4 # i gets shifted 2 bits to the right in every iteration 
                if secondmajorana == -1
                    if opi == 1 # so we got an X operator
                        secondmajorana = 2k # so the 2nd majoran op is of the 1st type
                        @debug "second operator is γ1 on site $k"
                    elseif opi == 2 # so we got an Y operator
                        secondmajorana = 2k+1 # so the 2nd majorana op is of the 2nd type
                        coeff *= -1
                        @debug "second operator is γ2 on site $k"
                    elseif opi == 3 # so we got a Z operator
                        firstmajorana = 2k+1    # so both majoranas act on the site
                        secondmajorana = 2k
                        @debug "second operator is γ2 on site $k"
                        @debug "first operator is γ1 on site $k"
                    else
                        @debug "nothing on site $k"
                    end
                elseif firstmajorana == -1
                    if opi == 1
                        firstmajorana = 2k+1
                        @debug "first operator is γ2 on site $k"
                    elseif opi == 2
                        firstmajorana = 2k
                        @debug "first operator is γ1 on site $k"
                    elseif opi == 3
                        coeff *= -1
                        @debug "nothing on site $k"
                    else
                        throw(NonFLOException("P contains terms that are not the product of two Majoranas"))
                    end
                elseif opi != 0
                    throw(NonFLOException("P contains terms that are not the product of two Majoranas"))
                else
                    @debug "nothing on site $k"
                end
                i = i >> 2
            end
            firstmajorana == -1 && throw(NonFLOException("P contains terms that are not the product of two Majoranas"))
            out[firstmajorana+1, secondmajorana+1] = -2 * real(coeff)
            out[secondmajorana+1, firstmajorana+1] =  2 * real(coeff)
        end
    end
    return out
end

# This function is not actually needed for the functionality of the simulator,
# but it is helpful for debugging purposes
"""
    majoranasquares2qubitbasis(H::AbstractMatrix)

Converts an ``2n×2n`` Majorana hamiltonian `H` into the full ``2^n×2^n`` hamiltonian
in the qubit basis.
"""
function majoranasquares2qubitbasis(H::AbstractMatrix)
    nq = size(H, 1) ÷ 2
    terms = [(1im * H[I] / 4) * (majoranaop(nq, I[1]) * majoranaop(nq, I[2]))
             for I in CartesianIndices(H)
             if abs(H[I]) >= √eps(abs(H[I]))]
    yaoblock = sum(terms)
    return mat(yaoblock)
end

"""
    qubit2majoranaevolution(U::AbstractMatrix, locs)

Turns a ``n`` qubit unitary `U` on the ``n`` qubits in `locs` into the corresponding
``SO(2n×2n)`` matrix for the evolution of the Majorana operators.
"""
function qubit2majoranaevolution(U::AbstractMatrix, locs)
    A = 1im * log(U)
    P = qubit2paulibasis(A)
    H = paulibasis2majoranasquares(P, locs)
    W = exp(H)
    return W
end

"""
    kron2majoranaindices(k::KronBlock)

Get the indices of majorana operators from a kronecker product of pauli operators
"""
function kron2majoranaindices(k::KronBlock)
    locs = collect(Iterators.flatten(k.locs))
    perm = sortperm(locs)
    areconsecutive(locs[perm]) || throw(NonFLOException("$k is not acting on consecutive qubits"))
    
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
                # @debug "second operator is γ1 on site $(locs[q])"
            elseif op == Y
                secondindex = 2locs[q]
                swap = !swap
                # @debug "second operator is γ2 on site $(locs[q])"
            elseif op == Z
                swap = !swap
            elseif op == I2
                # @debug "nothing on site $(locs[q])"
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

kron2majoranaindices(k::PutBlock{2,1,ZGate}) = (2k.locs[1]-1, 2k.locs[1])

"""
    yaoham2majoranasquares(::Type{T}=Float64, yaoham::AbstracBlock{2})

Convert a hamiltonian written as a YaoBlock into the corresponding 
``2n×2n`` majorana hamiltonian.
"""
function yaoham2majoranasquares(::Type{T}, yaoham::Add{2}) where {T<:Real}
    ham = zeros(T, 2nqubits(yaoham), 2nqubits(yaoham))
    for k in yaoham
        if k isa Scale
            i1, i2 = kron2majoranaindices(k.content)
            ham[i1,i2] += 2k.alpha
            ham[i2,i1] -= 2k.alpha
        elseif k isa KronBlock
            i1, i2 = kron2majoranaindices(k)
            ham[i1,i2] += 2
            ham[i2,i1] -= 2
        elseif k isa Add
            ham += yaoham2majoranasquares(T, k)
        elseif k isa PutBlock{2,1,ZGate} 
            i1, i2 = kron2majoranaindices(k)
            ham[i1,i2] += 2
            ham[i2,i1] -= 2
        else
            throw(NonFLOException("$k not recognized as a square of Majorana operators"))
        end
    end
    return ham
end

function yaoham2majoranasquares(::Type{T}, yaoham::KronBlock{2}) where {T<:Real}
    ham = spzeros(T, 2nqubits(yaoham), 2nqubits(yaoham))
    i1, i2 = kron2majoranaindices(yaoham)
    ham[i1,i2] = 2
    ham[i2,i1] = -2
    return ham
end

function yaoham2majoranasquares(::Type{T}, yaoham::PutBlock{2,1,ZGate}) where {T<:Real}
    ham = spzeros(T, 2nqubits(yaoham), 2nqubits(yaoham))
    i1, i2 = 2yaoham.locs[1]-1, 2yaoham.locs[1]
    ham[i1,i2] = 2
    ham[i2,i1] = -2
    return ham
end

function yaoham2majoranasquares(::Type{T}, yaoham::Scale) where {T<:Real}
    return yaoham.alpha * yaoham2majoranasquares(T, yaoham.content)
end

yaoham2majoranasquares(yaoham) = yaoham2majoranasquares(Float64, yaoham)

# -------------------------------------
# Utilities to generate random matrices
# -------------------------------------

"""
    random_unit_vector([Float64, ] n, N=n)

Generate a uniformly random vector on the `n`-sphere embedded in ℝ^N and 
return it and the sign of its first entry.
"""
function random_unit_vector(::Type{T}, n, N=n) where {T}
    v = [randn(T, n); zeros(N-n)]
    n = norm(v)
    s = sign(v[1])
    v[1] += s * n
    return v ./ norm(v), s
end

random_unit_vector(n) = random_unit_vector(Float64, n)


"""
    random_orthogonal_matrix([Float64,] n)

Generate a Haar random matrix in ``O(n)`` with element type `Float64` 
using the algorithm described in https://arxiv.org/pdf/math-ph/0609050.pdf
"""
function random_orthogonal_matrix(::Type{T}, n) where {T}
    out = Matrix{T}(I(n))
    for i in 2:n 
        v, s = random_unit_vector(T, i, n)
        out .*= -s
        out .-= 2v .* (v' * out)
    end
    rand(Bool) && (out[:,1] .*= -1) # randomize the determinant
    return out
end

random_orthogonal_matrix(n) = random_orthogonal_matrix(Float64, n)
