using CUDA
using CUDA.CUSPARSE
using KernelAbstractions
using ExponentialUtilities

export cu, cpu

const CuMajoranaReg{T} = MajoranaReg{T,CuArray{T,2,CUDA.DeviceMemory}} where {T}
const CPUMajoranaReg{T} = MajoranaReg{T,Matrix{T}} where {T}

CUDA.cu(reg::CPUMajoranaReg{T}) where T = MajoranaReg(reg.state |> CUDA.cu)
Yao.cpu(reg::CuMajoranaReg{T}) where T = MajoranaReg(Array(reg.state))


# Fallback for generic matrix gates. This is _much_ slower than the explicitely
# defined previous gates, so you do profit from defining unsafe_apply! and 
# instruct! for your fancy new FLO gate like I did before.
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


# Time-Evolution block specialisations
# ------------------------------------
function YaoBlocks.unsafe_apply!(reg::CuMajoranaReg, b::TimeEvolution)
    H = yaoham2majoranasquares(b.H) |> cu
    H .*= b.dt
    expH = exponential!(H) # maybe ExponetialUtilities can speed this up on GPU?
    reg.state .= expH * reg.state
    return reg
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
    kernel!(M, i, pi, ni, ndrange=(size(M,1)-2i, size(M,1)-2i-1))
    return M
end


# Expect specialisations
# ----------------------
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

# Autodiff support
# ----------------
function Yao.AD.backward_params!(st::Tuple{<:CuMajoranaReg,<:CuMajoranaReg},
                                 block::Rotor, collector)
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

function Yao.AD.expect_g(op::AbstractAdd, in::MajoranaReg)
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




# Benchmark functions
# -------------------
function cpu_benchmark_fun(n)
    reg = rand(2n, 2n) |> MajoranaReg
    covmat = covariance_matrix(reg)
end

function gpu_benchmark_fun(n)
    reg = CUDA.rand(2n, 2n) |> MajoranaReg
    covmat = CUDA.@sync covariance_matrix(reg)
end

function cpu_benchmark_fun(n)
    reg = rand(2n, 2n) |> MajoranaReg
    i = rand(1:n)
    pi = rand()
    ni = rand(Bool)
    #covmat = covariance_matrix(reg)
    update_covariance_matrix!(reg.state, i, pi, ni)
end

function gpu_benchmark_fun(n)
    reg = CUDA.rand(2n, 2n) |> MajoranaReg
    i = rand(1:n)
    pi = rand()
    ni = rand(Bool)
    #covmat = CUDA.@sync covariance_matrix(reg)
    CUDA.@sync FLOYao.update_covariance_matrix!(reg.state, i, pi, ni)
end
#=
function update_covariance_matrix!(M::CuMatrix, i, pi, ni)
    n = size(M,1)
    @inline function kernel!(M)
        p, q = blockIdx().x + 2i, threadIdx().x + 2i + 1
        offset = (-1)^ni * (M[2i-1,q] * M[2i,p] - M[2i-1,p] * M[2i,q]) / (2pi)
        M[p,q] += ifelse(q>p, offset, zero(eltype(M)))
        return nothing
    end
    @cuda threads=n-2i-1 blocks=n-2i kernel!(M)
    return M
end
=#


