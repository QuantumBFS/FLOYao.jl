"""
Support for fermionic linear optics circuits in Yao.jl. This is broadly based
on [Classical simulation of noninteracting-fermion quantum circuits](https://arxiv.org/abs/quant-ph/0108010),
but everything you need to know for the exact mathematical details of this
implementation is also found in my wiki in `simulating-fermionic-linear-optic-circuits`.
"""
module FLOYao

using LinearAlgebra
using Yao
using SparseArrays
using Random

export MajoranaReg

include("utils.jl")
include("majorana_reg.jl")
include("instruct.jl")
include("apply_composite.jl")
include("expect.jl")
include("measure.jl")
include("auto_diff.jl")

end # module


