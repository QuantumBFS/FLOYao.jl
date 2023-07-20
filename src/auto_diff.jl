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

const Rotor{T} = Union{RotationGate{2,T},
                       PutBlock{2,<:Any,<:RotationGate{<:Any,T}},
                       PutBlock{2,<:Any,PauliGate}
                      }

function Yao.AD.expect_g(op::AbstractAdd, in::MajoranaReg)
    ham = yaoham2majoranasquares(op)
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

function Yao.AD.backward_params!(st::Tuple{<:MajoranaReg,<:MajoranaReg},
                                 block::Rotor, collector)
    out, outδ = st
    ham = Yao.AD.generator(block)
    majoranaham = yaoham2majoranasquares(ham)
    g = outδ.state ⋅ (majoranaham * out.state) / 4
    pushfirst!(collector, g)
    return nothing
end

function Yao.AD.backward_params!(st::Tuple{<:MajoranaReg,<:MajoranaReg},
                                 block::TimeEvolution, collector)
    out, outδ = st
    ham = block.H
    majoranaham = yaoham2majoranasquares(ham)
    g = outδ.state ⋅ (majoranaham * out.state) / 2
    pushfirst!(collector, g)
    return nothing
end

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

function Yao.AD.expect_g(op::AbstractAdd, circuit::Pair{<:MajoranaReg,<:AbstractBlock})
    reg, c = circuit
    out = copy(reg) |> c
    outδ = Yao.AD.expect_g(op, out)
    paramsδ = [] 
    in, inδ = Yao.AD.apply_back!((out, outδ), c, paramsδ)
    # multiplying by 1 is a weird trick to  convert from Vector{Any} to Vector{Float}
    return inδ => 1paramsδ
end

