# Adding support for custom gates

Natively, the only FLO gates that come already shipped with `FLOYao.jl` are these
[Supported gates](@ref). But there are many more FLO gates,
one being for example the `FSWAP` gate which swaps to qubits while making sure
to preserve the fermionic commutation relations

```julia
@const_gate FSWAP::ComplexF64 = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 -1]
```

If a gate defines a matrix representation, as we just did for the `FSWAP`gate,
`FLOYao` supports them out of the box by manually checking if they are a FLO
gate and then computing its matrix representation in the Majorana basis. But
this method is fairly slow---though still poly-time and memory---compared to
directly implementing `unsafe_apply!(::MajoranaReg, ::YourBlock)` and
`instruct!(::MajoranaReg, ::YourBlock)` and will warn you accordingly:


```julia
nq = 4
fswap = put(nq, (1, 2) => FSWAP)
mreg = FLOYao.zero_state(nq)
mreg |> put(nq, 2 => X)
mreg |> fswap
```

    ┌ Warning: Calling manual instruct!(MajoranaReg{Float64}(4), ComplexF64[1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im -1.0 + 0.0im], (1, 2)).
    │ You can greatly speed up your FLO gates by exactly implementing unsafe_apply!()
    │ and instruct!() for them. See FLOYao/src/instruct.jl and  FLOYao/src/apply_composite.jl
    │ for how to do that.
    └ @ FLOYao ~/.julia/.../FLOYao/src/instruct.jl:56

    MajoranaReg{Float64} with 4 qubits:
     -2.35415e-16  -4.12493e-16  -1.0          …   0.0   0.0   0.0   0.0
      2.46746e-16  -5.5708e-16   -1.26504e-16      0.0   0.0   0.0   0.0
     -1.0          -1.17708e-16   2.55988e-16      0.0   0.0   0.0   0.0
     -1.85286e-16  -1.0           2.44068e-16      0.0   0.0   0.0   0.0
     -0.0          -0.0          -0.0             -1.0  -0.0  -0.0  -0.0
     -0.0          -0.0          -0.0          …  -0.0  -1.0  -0.0  -0.0
     -0.0          -0.0          -0.0             -0.0  -0.0  -1.0  -0.0
     -0.0          -0.0          -0.0             -0.0  -0.0  -0.0  -1.0


Now, before we fix these warnings, let's see how long the current implementation takes:


```julia
using BenchmarkTools
using Suppressor # we don't want to get all the warnings when benchmarking
@benchmark @suppress apply!($mreg, $fswap)
```

    BenchmarkTools.Trial: 6524 samples with 1 evaluation.
     Range (min … max):  707.996 μs …   3.773 ms  ┊ GC (min … max): 0.00% … 73.67%
     Time  (median):     727.015 μs               ┊ GC (median):    0.00%
     Time  (mean ± σ):   761.750 μs ± 225.610 μs  ┊ GC (mean ± σ):  2.31% ±  6.25%
    
      ▆█▇▆▅▄▄▄▃▃▂▂▁▁                                                ▂
      ███████████████▇██▇▇▃▅▅▄▅▅▄▃▄▁▃▃▁▁▃▃▃▃▁▃▁▃▁▁▁▃▁▄▁▁▁▃▃▆▆▆▆▆▆▄▄ █
      708 μs        Histogram: log(frequency) by time        1.2 ms <
    
     Memory estimate: 338.53 KiB, allocs estimate: 495.


To find out what the matrix representation of the `FSWAP` gate in the Majorana
basis is, it is easiest to retrace what is happening inside
`instruct!(::MajoranaReg, ::AbstractMatrix, locs)`. You can use


```julia
@which instruct!(mreg, mat(FSWAP), (1,2))
```

    instruct!(reg::MajoranaReg, gate::AbstractMatrix, locs) in FLOYao at ~/.julia/.../FLOYao/src/instruct.jl:49


to find the location of the corresponding code. Now let's copy-paste what we found there:

```julia
W = FLOYao.qubit2majoranaevolution(Matrix(fswap.content), fswap.locs)
```

    4×4 Matrix{Float64}:
     -2.35415e-16  -4.12493e-16  -1.0           0.0
      2.46746e-16  -5.5708e-16   -1.26504e-16  -1.0
     -1.0          -1.17708e-16   2.55988e-16  -2.38988e-16
     -1.85286e-16  -1.0           2.44068e-16   2.43374e-16


```julia
matlocs = 2*(fswap.locs[1]-1)+1:2(fswap.locs[end])
```
    1:4

this matrix gets left-multiplied to the columns `1:4` in the last line of
`FLOYao.majorana_unsafe_apply!(::MajoranaReg, ::PutBlock)`. So we can instead
implement the action of an `FSWAP` gate on a `MajoranaReg` directly as follows:

```julia
function YaoBlocks.unsafe_apply!(reg::MajoranaReg, b::PutBlock{2,2,FSWAPGate})
    FLOYao.areconsecutive(b.locs) || throw(NonFLOException("FSWAP must act on consecutive qubits"))
    instruct!(reg, Val(:FSWAP), b.locs)
end

function Yao.instruct!(reg::MajoranaReg, ::Val{:FSWAP}, locs::Tuple)
    i1, i2 = locs
    row1, row2 = reg.state[2i1-1,:], reg.state[2i1,:]
    row3, row4 = reg.state[2i2-1,:], reg.state[2i2,:]
    reg.state[2i1-1,:] .=  .-row3
    reg.state[2i1,:] .=  .-row4
    reg.state[2i2-1,:] .=  .-row1
    reg.state[2i2,:] .=  .-row2
    return reg
end
```

```julia
@benchmark apply!($mreg, $fswap)
```

    BenchmarkTools.Trial: 10000 samples with 676 evaluations.
     Range (min … max):  183.416 ns …   4.023 μs  ┊ GC (min … max): 0.00% … 93.94%
     Time  (median):     198.760 ns               ┊ GC (median):    0.00%
     Time  (mean ± σ):   224.080 ns ± 242.728 ns  ┊ GC (mean ± σ):  9.29% ±  8.02%
    
      ▇ ▆ █▁▅ ▁                                                      
      █▅███████▆▆▅▄▃▃▃▂▃▂▂▂▂▂▂▂▂▂▂▂▁▂▂▂▂▂▂▂▁▁▁▂▂▁▁▁▂▂▂▁▂▁▂▂▂▂▂▂▂▂▂▂ ▃
      183 ns           Histogram: frequency by time          362 ns <
    
     Memory estimate: 512 bytes, allocs estimate: 4.

Which is indeed a significant speed-up!
