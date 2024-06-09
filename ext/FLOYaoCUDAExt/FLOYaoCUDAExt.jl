module FLOYaoCUDAExt

using CUDA
using CUDA.CUSPARSE
using ExponentialUtilities
using FLOYao
using KernelAbstractions
using Yao
using LinearAlgebra

import CUDA: cu
import Yao: cpu

import FLOYao: majoranaop, cuproduct_state, cuzero_state, cuone_state
import FLOYao: one_state!, zero_state!, product_state!, majorana2arrayreg
import FLOYao: covariance_matrix, update_covariance_matrix!, sample,
    majorana_expect, bitstring_probability, yaoham2majoranasquares, fast_overlap


const CuMajoranaReg{T} = MajoranaReg{T,CuArray{T,2,CUDA.DeviceMemory}} where {T}
const CPUMajoranaReg{T} = MajoranaReg{T,Matrix{T}} where {T}

# Register creation and transfer
# ------------------------------
cu(reg::CPUMajoranaReg{T}) where T = MajoranaReg(reg.state |> cu)
cpu(reg::CuMajoranaReg{T}) where T = MajoranaReg(Array(reg.state))


# majorana_reg.jl extensions
# --------------------------
function cuzero_state(::Type{T}, n::Integer) where {T}
    state = CUDA.zeros(T, 2n, 2n)
    state[diagind(state)] .= one(T)
    return MajoranaReg(state)
end

cuzero_state(n) = cuzero_state(Float32, n)

@kernel function one_state_kernel!(state)
    i, j = @index(Global, NTuple)
    if i != j
        state[i,j] = 0
    else
        state[i,j] = (i % 4 == 1 ||  i % 4 == 0) ? 1 : -1
    end
end

function one_state!(reg::CuMajoranaReg{T}) where {T}
    backend = KernelAbstractions.get_backend(reg.state)
    kernel! = one_state_kernel!(backend)
    kernel!(reg.state, ndrange=size(reg.state))
    return reg
end

function cuone_state(::Type{T}, n::Integer) where {T}
    reg = CUDA.CuMatrix{T}(undef, 2n, 2n) |> MajoranaReg
    return one_state!(reg)
end

cuone_state(n) = cuone_state(Float32, n)

function cuproduct_state(::Type{T}, bit_str::DitStr{2,N,IT}) where {N,T,IT}
    reg = MajoranaReg(CuMatrix{T}(undef, 2N, 2N))
    product_state!(reg, bit_str)
    return reg
end

cuproduct_state(bit_str) = cuproduct_state(Float32, bit_str)

# vector input
function cuproduct_state(::Type{T}, bit_configs::AbstractVector) where {T}
    return cuproduct_state(T, DitStr{2}(bit_configs))
end

# integer input
function cuproduct_state(::Type{T}, nbits::Int, val::Integer) where {T}
    return cuproduct_state(T, DitStr{2,nbits}(val))
end

function majorana2arrayreg(reg::CuMajoranaReg)
    nq = nqubits(reg)

    # praying here, that the rand_state has non-zero overlap with the 
    # new vacuum state.
    # This first bit gets areg into the new vacuum by piping it through
    # the projector U|Ω⟩⟨Ω|U^† 
    areg = Yao.curand_state(Complex{eltype(reg)}, nq)
    for i in 1:nq
        γ_i1 = majoranaop(nq, reg.state[:,2i-1] |> Vector)
        γ_i2 = majoranaop(nq, reg.state[:,2i] |> Vector)
        circuit = 1im*chain(nq, γ_i2, γ_i1) + igate(nq)
        areg |> circuit
    end
    normalize!(areg)
    return areg
end


# instruct.jl extensions
# ----------------------
function Yao.instruct!(reg::CuMajoranaReg, gate::AbstractMatrix, locs)
    # unless the gate is a tensor product of Z's (but _not_ using Yaos kron
    # block a gate acting on non-consecutive qubits will never be a 
    # FLO gate
    areconsecutive(locs) || throw(NonFLOException("$gate on $locs is not a FLO gate"))
    W = qubit2majoranaevolution(gate, locs) |> cu
    matlocs = (2*(locs[1]-1)+1:2(locs[end]))
    @warn """Calling manual instruct!($reg, $gate, $locs).
    You can greatly speed up your FLO gates by exactly implementing unsafe_apply!()
    and instruct!() for them. See FLOYao/src/instruct.jl and  FLOYao/src/apply_composite.jl
    for how to do that.
    """
    reg.state[matlocs,:] .= W * reg.state[matlocs,:]
end


# apply_composite.jl extensions
# -----------------------------
function YaoBlocks.unsafe_apply!(reg::CuMajoranaReg, b::TimeEvolution)
    H = yaoham2majoranasquares(b.H) |> cu
    H .*= b.dt
    expH = exponential!(H)
    reg.state .= expH * reg.state
    return reg
end


# measure.jl extensions
# ---------------------
# Can probably make this even faster by writing out what
# happens in index notation and doing it _all_ in the kernel
@kernel function swap_and_sign_kernel!(A)
    i, j = @index(Global, NTuple)
    j = 2 * j
    A[i,j], A[i,j-1] = A[i,j-1], -A[i,j]
end

function covariance_matrix(reg::CuMajoranaReg)
    nq = nqubits(reg)
    state = reg.state |> copy
    backend = KernelAbstractions.get_backend(state)
    kernel! = swap_and_sign_kernel!(backend)
    kernel!(state, ndrange=(size(state, 1), size(state,2) ÷ 2))
    return state * reg.state'
end

@kernel function update_covariance_matrix_kernel!(M, i, pi, ni)
    p, q = @index(Global, NTuple)
    p += 2i
    q += 2i + 1
    offset = (-1)^ni * (M[2i-1,q] * M[2i,p] - M[2i-1,p] * M[2i,q]) / (2pi)
    M[p,q] += ifelse(q>p, offset, zero(eltype(M)))
end

function update_covariance_matrix!(M::CuMatrix, i, pi, ni)
    backend = KernelAbstractions.get_backend(M)
    kernel! = update_covariance_matrix_kernel!(backend)
    return if size(M, 1) <= 2i # nothing to do for last qubit
        M
    else
        kernel!(M, i, pi, ni, ndrange=(size(M,1)-2i, size(M,1)-2i-1))
    end
end


@inline function sample(covmat::CuMatrix, locs, ids,
                        rng=Random.GLOBAL_RNG)
    out = BigInt(0)
    nq = length(locs)
    for i in 1:nq
        pi = CUDA.@allowscalar (1 + covmat[2locs[ids[i]]-1,2locs[ids[i]]]) / 2
        ni = rand(rng) > pi
        out += ni * 2^(ids[i]-1)
        update_covariance_matrix!(covmat, i, ni ? 1-pi : pi, ni)
    end
    return out
end

# a slightly faster version, when all qubits get sampled in their normal
# order
@inline function sample(covmat::CuMatrix, rng=Random.GLOBAL_RNG)
    out = BigInt(0)
    nq = size(covmat,1) ÷ 2
    for i in 1:nq
        pi = CUDA.@allowscalar (1 + covmat[2i-1,2i]) / 2
        ni = rand(rng) > pi
        out += ni * 2^(i-1)
        update_covariance_matrix!(covmat, i, ni ? 1-pi : pi, ni)
    end
    return out
end

# expect.jl extensions
# --------------------
# For actual performance, it would be better if we could construct block on the GPU!
# Also, the actual logic in majorana_expect can probably be made faster on a GPU
function majorana_expect(block::AbstractMatrix, reg::CuMajoranaReg)
    block = block |> cu

    expval = sum(1:nqubits(reg), init=zero(eltype(reg.state))) do i
        @views reg.state[:,2i] ⋅ (block * reg.state[:,2i-1]) - reg.state[:,2i-1] ⋅ (block * reg.state[:,2i])
    end

    offset_expval = sum(1:nqubits(reg), init=zero(eltype(reg.state))) do i
        @views - reg.state[:,2i] ⋅ (block * reg.state[:,2i]) - reg.state[:,2i-1] ⋅ (block * reg.state[:,2i-i])
    end |> imag

    return (- expval - offset_expval) / 4
end

function bitstring_probability(reg::CuMajoranaReg{T}, bit_string::DitStr{2,N,ST}) where {T,N,ST}
    @assert nqubits(reg) == N
    p = one(T)
    M = covariance_matrix(reg)
    for i in 1:N
        ni = bit_string[i]
        pi = CUDA.@allowscalar (1 + (-1)^ni * M[2i-1,2i]) / 2
        pi ≈ 0 && return zero(T)
        p *= pi
        update_covariance_matrix!(M, i, pi, ni)
    end
    return p > zero(T) ? p : zero(T) # floating point issues can cause very small probabilities to get negative.
end

# auto_diff.jl extensions
# -----------------------
function Yao.AD.backward_params!(st::Tuple{<:CuMajoranaReg,<:CuMajoranaReg},
                                 block::FLOYao.Rotor, collector)
    out, outδ = st
    ham = Yao.AD.generator(block)
    majoranaham = yaoham2majoranasquares(ham) |> cu
    g = fast_overlap(outδ.state, majoranaham, out.state) / 4
    pushfirst!(collector, g)
    return nothing
end

function Yao.AD.backward_params!(st::Tuple{<:CuMajoranaReg,<:CuMajoranaReg},
                                 block::TimeEvolution, collector)
    out, outδ = st
    ham = block.H
    majoranaham = yaoham2majoranasquares(ham) |> cu
    g = fast_overlap(outδ.state, majoranaham, out.state) / 2
    pushfirst!(collector, g)
    return nothing
end

function Yao.AD.expect_g(op::AbstractAdd, in::CuMajoranaReg)
    ham = yaoham2majoranasquares(op) |> cu
    inδ = copy(in)
    inδ.state .= ham * inδ.state
    for i in 1:nqudits(in)
        ψ1, ψ2 = inδ.state[:,2i-1], inδ.state[:,2i]
        inδ.state[:,2i-1] .= -ψ2
        inδ.state[:,2i] .= ψ1
    end
    inδ.state[:,1:end] .*= -1
    return inδ
end

end
