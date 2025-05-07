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

# -------------------------------------
# Getting derivatives w.r.t. parameters
# -------------------------------------
const Rotor{T} = Union{RotationGate{2,T},
                       PutBlock{2,<:Any,<:RotationGate{<:Any,T}},
                       PutBlock{2,<:Any,PauliGate}
                      }

function Yao.AD.backward_params!(st::Tuple{<:MajoranaReg,<:MajoranaReg},
                                 block::Rotor, collector)
    out, outδ = st
    ham = Yao.AD.generator(block)
    majoranaham = yaoham2majoranasquares(ham)
    g = fast_overlap(outδ.state, majoranaham, out.state) / 4
    pushfirst!(collector, g)
    return nothing
end

function Yao.AD.backward_params!(st::Tuple{<:MajoranaReg,<:MajoranaReg},
                                 block::TimeEvolution, collector)
    out, outδ = st
    ham = block.H
    majoranaham = yaoham2majoranasquares(ham)
    g = fast_overlap(outδ.state, majoranaham, out.state) / 2
    pushfirst!(collector, g)
    return nothing
end

# --------------------------
# Reverse propagating states
# --------------------------
function Yao.AD.apply_back(st::Tuple{<:MajoranaReg,<:MajoranaReg}, block::AbstractBlock)
    paramsδ = []
    in, inδ = Yao.AD.apply_back!(st, block, paramsδ)
    return (in, inδ), paramsδ
end

function Yao.AD.apply_back!(st::Tuple{<:MajoranaReg,<:MajoranaReg},
                            block::Union{PutBlock{2,<:Any,BT},TimeEvolution{2,<:Any,BT}},
                            collector) where {BT}
    out, outδ = st
    adjblock = block'
    if nparameters(block) != 0
        Yao.AD.backward_params!((out, outδ), block, collector)
    end
    in = apply!(out, adjblock)
    inδ = apply!(outδ, adjblock)
    return (in, inδ)
end

# -------------------------------
# Gradients of expectation values
# -------------------------------
function Yao.AD.expect_g(op::AbstractAdd, in::MajoranaReg)
    ham = yaoham2majoranasquares(op)
    inδ = copy(in)
    inδ.state .= ham * inδ.state
    @inbounds for i in 1:nqudits(in)
        for k in 1:size(inδ.state, 1)
            inδ.state[k,2i-1], inδ.state[k,2i] = -inδ.state[k,2i], inδ.state[k,2i-1]
        end
    end
    inδ.state[:,1:end] .*= -1
    return inδ
end


function Yao.AD.expect_g(op::AbstractAdd, circuit::Pair{<:MajoranaReg,<:AbstractBlock})
    reg, c = circuit
    out = copy(reg) |> c
    outδ = Yao.AD.expect_g(op, out)
    paramsδ = [] 
    in, inδ = Yao.AD.apply_back!((out, outδ), c, paramsδ)
    # multiplying by 1 is a weird trick to  convert from Vector{Any} to Vector{Float}
    return inδ => 1paramsδ
end

# -----------------------
# Gradients of fidelities
# -----------------------
function (::Adjoint{Any,typeof(Yao.fidelity)})(
    reg1::Union{MajoranaReg,Pair{<:MajoranaReg,<:AbstractBlock}},
    reg2::Union{MajoranaReg,Pair{<:MajoranaReg,<:AbstractBlock}},
)
    Yao.AD.fidelity_g(reg1, reg2)
end

function Yao.AD.fidelity_g(
        reg1::Union{MajoranaReg,Pair{<:MajoranaReg,<:AbstractBlock}},
        reg2::Union{MajoranaReg,Pair{<:MajoranaReg,<:AbstractBlock}}
    )

    if reg1 isa Pair
        in1, c1 = reg1
        out1 = copy(in1) |> c1
    else
        out1 = reg1
    end

    if reg2 isa Pair
        in2, c2 = reg2
        out2 = copy(in2) |> c2
    else
        out2 = reg2
    end

    R = out1.state' * out2.state
    # different parity => early out
    if det(R) < 0
        if reg1 isa Pair
            res1 = MajoranaReg(zero(R)) => zero(parameters(c1))
        else
            res1 = MajoranaReg(zero(R))
        end

        if reg2 isa Pair
            res2 = MajoranaReg(zero(R)) => zero(parameters(c2))
        else
            res2 = MajoranaReg(zero(R))
        end

        return (res1, res2)
    end

    fidelityδ = fidelity_gradient(R)
    out1δ = MajoranaReg(out2.state * fidelityδ')
    out2δ = MajoranaReg(out1.state * fidelityδ)

    if reg1 isa Pair
        (_, in1δ), params1δ = Yao.AD.apply_back((out1, out1δ), c1)
        res1 = in1δ => 1params1δ
    else
        res1 = out1δ
    end

    if reg2 isa Pair
        (_, in2δ), params2δ = Yao.AD.apply_back((out2, out2δ), c2)
        res2 = in2δ => 1params2δ
    else
       res2 = out2δ
    end

    return (res1, res2)
end
