# PutBlock specialisations
# ------------------------
# The absolute fallback. Probably not a good idea
function YaoBlocks.unsafe_apply!(reg::MajoranaReg, pb::PutBlock)
    instruct!(reg, mat(pb.content), pb.locs)
    return reg
end

for G in [:X, :Y, :Z, :T, :Tdag]
    GT = Expr(:(.), :ConstGate, QuoteNode(Symbol(G, :Gate)))
    @eval function Yao.unsafe_apply!(reg::MajoranaReg, pb::PutBlock{2,1,<:$GT})
        instruct!(reg, Val($(QuoteNode(G))), pb.locs)
        return reg
    end
end

for (G, GT) in [(:Rz, :(RotationGate{2,T,ZGate})),
                (:Shift, :(ShiftGate{T}))]
    @eval function Yao.unsafe_apply!(reg::MajoranaReg, pb::PutBlock{2,1,$GT}) where {T}
        instruct!(reg, Val($(QuoteNode(G))), pb.locs, pb.content.theta)
        return reg
    end
end

for GT in [:(RotationGate{2,T,XGate}),
           :(RotationGate{2,T,YGate})]
    @eval function Yao.unsafe_apply!(reg::MajoranaReg, pb::PutBlock{2,1,$GT}) where {T}
        throw(NonFLOException("Rx and Ry are not a FLO gates"))
        return reg
    end
end

# KronBlock specialisations
# -------------------------
# The default solution with fallback
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

function YaoBlocks.unsafe_apply!(reg::MajoranaReg, k::KronBlock)
    backup = copy(reg.state)
    try   # first try just applying everything block by block
        nq = nqubits(k)
        for (locs, block) in zip(k.locs, k.blocks)
            YaoBlocks.unsafe_apply!(reg, put(nq, Tuple(locs) => block))
        end
    catch # and if that fails see if maybe the whole product becomes a FLO gate
        reg.state .= backup
        # locs = sort!(collect(Iterators.flatten(k.locs)))
        locs = collect(Iterators.flatten(k.locs))
        instruct!(reg, reducedmat(k.content), locs)
    end
    return reg
end

const PauliGate = Union{I2Gate,XGate,YGate,ZGate}
const PauliKronBlock = KronBlock{2,N,<:NTuple{N,PauliGate}} where {N}
function YaoBlocks.unsafe_apply!(reg::MajoranaReg, k::PauliKronBlock)
    i1, i2 = kron2majoranaindices(k)
    reg.state[i1,:] .*= -1
    reg.state[i2,:] .*= -1
    return reg
end

const RGate = RotationGate{2,<:Real,<:PauliKronBlock}
function YaoBlocks.unsafe_apply!(reg::MajoranaReg, b::RGate)
    i1, i2 = kron2majoranaindices(b.block)
    s, c = sincos(b.theta)
    ψ1, ψ2 = reg.state[i1,:], reg.state[i2,:]
    reg.state[i1,:] .= c .* ψ1 .+ s .* ψ2
    reg.state[i2,:] .= c .* ψ2 .- s .* ψ1
    return reg
end

# Repeat Block specialisations
# ----------------------------
function YaoBlocks.unsafe_apply!(reg::MajoranaReg, rp::RepeatedBlock)
    nq = nqubits(reg)
    for addr in rp.locs
        YaoBlocks.unsafe_apply!(reg, put(nq, Tuple(addr:addr+nqudits(rp.content)-1) =>  rp.content))
    end
    return reg
end

for G in [:X, :Y, :Z, :T, :Tdag]
    GT = Expr(:(.), :ConstGate, QuoteNode(Symbol(G, :Gate)))
    @eval function Yao.unsafe_apply!(reg::MajoranaReg, rp::RepeatedBlock{2,N,<:$GT}) where {N}
        for addr in rp.locs
            instruct!(reg, Val($(QuoteNode(G))), Tuple(addr))
        end
        return reg
    end
end
