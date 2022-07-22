# ------------------------------
# Applying gates to the register
# ------------------------------

function Yao.instruct!(reg::MajoranaReg, ::Val{:Z}, locs::Tuple)
    loc = locs[1]
    reg.state[2loc-1:2loc,:] .*= -1
    return reg
end

function Yao.instruct!(reg::MajoranaReg, ::Val{:X}, locs::Tuple)
    loc = locs[1]
    reg.state[2loc:end,:] .*= -1
    return reg
end

function Yao.instruct!(reg::MajoranaReg, ::Val{:Y}, locs::Tuple)
    loc = locs[1]
    reg.state[2loc+1:end,:] .*= -1
    reg.state[2loc-1,:] .*= -1
    return reg
end

function Yao.instruct!(reg::MajoranaReg, ::Val{:Rz}, locs::Tuple, theta)
    loc = locs[1]
    c, s = cos(theta), sin(theta)
    ψ1, ψ2 = reg.state[2loc-1,:], reg.state[2loc,:]
    reg.state[2loc-1,:] .= c .* ψ1 .+ s .* ψ2
    reg.state[2loc,:] .= c .* ψ2 .- s .* ψ1
    return reg
end

const PauliGate = Union{I2Gate,XGate,YGate,ZGate}
const PauliKronBlock = KronBlock{2,N,<:NTuple{N,PauliGate}} where {N}
function YaoBlocks.unsafe_apply!(reg::MajoranaReg, k::PauliKronBlock) where {N}
    i1, i2 = kron2majoranaindices(k)
    reg.state[i1,:] .*= -1
    reg.state[i2,:] .*= -1
    return reg
end

const RGate = RotationGate{2,<:Real,<:PauliKronBlock}
function YaoBlocks.unsafe_apply!(reg::MajoranaReg, b::RGate)
    i1, i2 = kron2majoranaindices(b.block)
    c, s = cos(b.theta), sin(b.theta)
    ψ1, ψ2 = reg.state[i1,:], reg.state[i2,:]
    reg.state[i1,:] .= c .* ψ1 .+ s .* ψ2
    reg.state[i2,:] .= c .* ψ2 .- s .* ψ1
    return reg
end

# Fallback for generic matrix gates. This is _much_ slower than the explicitely
# defined previous gates, so you do profit from defining unsafe_apply! and 
# instruct! for your fancy new FLO gate like I did before.
function Yao.instruct!(reg::MajoranaReg, gate::AbstractMatrix, locs)
    # unless the gate is a tensor product of Z's (but _not_ using Yaos kron
    # block a gate acting on non-consecutive qubits will never be a 
    # FLO gate
    areconsecutive(locs) || throw(NonFLOException("$gate on $locs is not a FLO gate"))
    W = qubit2majoranaevolution(gate, locs)
    matlocs = (2*(locs[1]-1)+1:2(locs[end]))
    @warn """Calling manual instruct!($reg, $gate, $locs).
    You can greatly speed up your FLO gates by exactly implementing unsafe_apply!()
    and instruct!() for them. See FLOYao/src/instruct.jl for how to do that.
    """
    reg.state[matlocs,:] .= W * reg.state[matlocs,:]
end

"""
    reducedmat(::Type{T}, k::KronBlock)

Matrix representation of `k` with only the qubits `k` acts on non-trivially
in the tensor product. (As opposed to Yao.mat(k) which gives the full 2^n×2^n
matrix)
"""
function reducedmat(::Type{T}, k::KronBlock{D,M}) where {T,D,M}
    M == 0 && return IMatrix{D^k.n,T}()
    return reduce(
        Iterators.reverse(subblocks(k)),
        init = [one(T)],
    ) do x, y
        kron(x, mat(T, y))
    end |> Matrix
end

reducedmat(b::KronBlock) = reducedmat(promote_type(ComplexF64, parameters_eltype(b)), b)

function majorana_unsafe_apply!(reg::MajoranaReg, pb::PutBlock)
    try 
        instruct!(reg, pb.content, pb.locs)
    catch e
        e isa NonFLOException || e isa MethodError || rethrow(e)
        W = qubit2majoranaevolution(Matrix(pb.content), pb.locs)
        matlocs = (2*(pb.locs[1]-1)+1:2(pb.locs[end]))
        reg.state[matlocs,:] .= W * reg.state[matlocs,:]
    end
end

function YaoBlocks.unsafe_apply!(reg::MajoranaReg, k::KronBlock)
    backup = copy(reg.state)
    try   # first try just applying everything block by block
        for (locs, block) in zip(k.locs, k.blocks)
            instruct!(reg, block, Tuple(locs))
        end
    catch # and if that fails see if maybe the whole product becomes a FLO gate
        reg.state .= backup
        # locs = sort!(collect(Iterators.flatten(k.locs)))
        locs = collect(Iterators.flatten(k.locs))
        W = qubit2majoranaevolution(reducedmat(k), locs)
        matlocs = (2*(locs[1]-1)+1:2(locs[end]))
        reg.state[matlocs,:] .= W * reg.state[matlocs,:] 
    end
    return reg
end

