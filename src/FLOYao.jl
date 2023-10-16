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

"""
A [Yao.jl](https://github.com/QuantumBFS/Yao.jl) backend to efficiently simulate
fermionic linear optics (FLO) circuits  based on
[Classical simulation of noninteracting-fermion quantum circuits](https://arxiv.org/abs/quant-ph/0108010)
and [Disorder-assisted error correction in Majorana chains](https://arxiv.org/abs/1108.3845).
FLO circuits are a class of quantum circuits that are closely related to
non-interacting fermions and can be efficiently simulated on classical
computers, similar to the way Clifford circuits can be efficiently classically
simulated, as is done in [YaoClifford.jl](https://github.com/QuantumBFS/YaoClifford.jl).

The goal of `FLOYao.jl` is that if you have code written in `Yao.jl` that only 
uses [FLO gates](https://quantumbfs.github.io/FLOYao.jl/stable/supported_gates/)
and other primitives that are efficiently simulatable in polynomial time and 
space, that you can simply replace your `AbstractArrayReg` with a `MajoranaReg`
and run exactly the same simulation, with the same code but exponentially faster.

A brief introduction to fermionic linear optics circuits is found in the 
[Documentation](https://yaoquantum.org/FLOYao.jl/stable/background/) and a more
in-depth introduction in e.g. the two papers linked above.
"""
module FLOYao

using LinearAlgebra
using Random
using SparseArrays
using Yao

export MajoranaReg

const PauliGate = Union{I2Gate,XGate,YGate,ZGate}
const PauliKronBlock = KronBlock{2,N,<:NTuple{N,PauliGate}} where {N}
const RGate = RotationGate{2,<:Real,<:PauliKronBlock}

include("utils.jl")
include("majorana_reg.jl")
include("instruct.jl")
include("apply_composite.jl")
include("expect.jl")
include("measure.jl")
include("auto_diff.jl")

end # module


