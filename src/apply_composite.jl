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
        instruct!(reg, reducedmat(k), locs)
    end
    return reg
end

# Defined in FLOYao.jl
# const PauliGate = Union{I2Gate,XGate,YGate,ZGate}
# const PauliKronBlock = KronBlock{2,N,<:NTuple{N,PauliGate}} where {N}
function YaoBlocks.unsafe_apply!(reg::MajoranaReg, k::PauliKronBlock)
    i1, i2 = kron2majoranasquare(k)
    reg.state[i1,:] .*= -1
    reg.state[i2,:] .*= -1
    return reg
end

# Defined in FLOYao.jl
# const RGate = RotationGate{2,<:Real,<:PauliKronBlock}
function YaoBlocks.unsafe_apply!(reg::MajoranaReg, rgate::RGate)
    i1, i2 = kron2majoranasquare(rgate.block)
    s, c = sincos(rgate.theta)
    for k in 1:size(reg.state, 2)
        ψ1, ψ2 = reg.state[i1,k], reg.state[i2,k]
        reg.state[i1,k] = c * ψ1 + s * ψ2
        reg.state[i2,k] = c * ψ2 - s * ψ1
    end
    return reg
end

function YaoBlocks.unsafe_apply!(reg::MajoranaReg, rpb::PutBlock{2,N,<:RGate}) where {N}
    areconsecutive(rpb.locs) || throw(NonFLOException("$(rpb.blocks) on $(rpb.locs) is not a FLO gate"))
    # goddamnit, 1-based indexing
    i1, i2 = 2 * (minimum(rpb.locs) - 1) .+ kron2majoranasquare(rpb.content.block)
    s, c = sincos(rpb.content.theta)
    for k in 1:size(reg.state, 2)
        ψ1, ψ2 = reg.state[i1, k], reg.state[i2, k]
        reg.state[i1, k] = c * ψ1 + s * ψ2
        reg.state[i2, k] = c * ψ2 - s * ψ1
    end
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

# Time-Evolution block specialisations
# ------------------------------------
function YaoBlocks.unsafe_apply!(reg::MajoranaReg, b::TimeEvolution)
    H = yaoham2majoranasquares(b.H)
    reg.state .= exp(b.dt .* H) * reg.state
    return reg
end
