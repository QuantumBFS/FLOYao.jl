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
Support for fermionic linear optics circuits in Yao.jl. This is broadly based
on [Classical simulation of noninteracting-fermion quantum circuits](https://arxiv.org/abs/quant-ph/0108010),
but everything you need to know for the exact mathematical details of this
implementation is also found in my wiki in `simulating-fermionic-linear-optic-circuits`.
"""
module FLOYao

using LinearAlgebra
using Random
using SparseArrays
using Yao

export MajoranaReg

include("utils.jl")
include("majorana_reg.jl")
include("instruct.jl")
include("apply_composite.jl")
include("expect.jl")
include("measure.jl")
include("auto_diff.jl")

end # module


