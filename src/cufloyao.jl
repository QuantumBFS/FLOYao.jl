using CUDA
using CUDA.CUSPARSE

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
#=
"Majorana op on `n` Majorana sites for orbital `ψ`"
function majoranaop(n, ψ::CuArray)
    ψ = ψ |> Vector
    return sum([ψ[i] * majoranaop(n, i) for i in eachindex(ψ)])
end
=#
# Can probably make this even faster by writing out what
# happens in index notation and doing it _all_ in the kernel
function covariance_matrix(reg::CuMajoranaReg)
    nq = nqubits(reg)
    state = reg.state |> copy
    threads = nq
    blocks = 2nq
    @inline function kernel(state)
        i, j = blockIdx().x, 2 * threadIdx().x
        tmp = state[i,j]
        state[i,j] = state[i,j-1]
        state[i,j-1] = -tmp
        nothing
    end
    @cuda threads=threads blocks=blocks kernel(state)
    return state * reg.state'
end

# This seems to be not really making fully use of the GPU and is barely faster
# than on the CPU. Also for _very_ large systems doing thousands of blocks 
# seems to be disallowed.
#=
    n = size(M,1)
    for p in 2i+1:n
        for q in p+1:n
            M[p,q] += (-1)^ni * M[2i-1,q] * M[2i,p] / (2pi)
            M[p,q] -= (-1)^ni * M[2i-1,p] * M[2i,q] / (2pi)
        end
    end
=#
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
