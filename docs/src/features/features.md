# Features

The core part of `FLOYao` is the a new register type, the `MajoranaReg`:

```@docs
MajoranaReg
```

## State initialization

There are several functions provided, to create [`MajoranaReg`](@ref)'s 
in different states. They are not exported from `FLOYao.jl` in order to avoid
name collisions with `Yao.jl`

```@docs
FLOYao.zero_state
FLOYao.zero_state_like
FLOYao.one_state
FLOYao.product_state
FLOYao.rand_state
```

as well as some functions to reset [`MajoranaReg`](@ref)'s to fixed states:

```@docs
FLOYao.zero_state!
FLOYao.one_state!
FLOYao.product_state!
```

## Applying gates
Application of [Supported gates](@ref) to [`MajoranaReg`](@ref)'s works 
as it does in [`Yao.jl`](https://docs.yaoquantum.org/stable/man/blocks.html#Base.:|%3E-Tuple{AbstractRegister,%20AbstractBlock})

```julia
reg = FLOYao.zero_state(2)
gate = put(2, 1 => X)

# using apply!
apply!(reg, gate)

# or using the pipe syntax
reg |> gate

# or the non-mutating version
result = apply(reg, gate)
```

## Measuring expectation values
The same goes measuring expectation values 
```@docs
expect
```
where (for now) `op` is restricted to be an observable quadratic in the 
Majorana operators (see [Known restrictions](@ref)). 
See the section about [`expect`](https://docs.yaoquantum.org/stable/man/blocks.html#YaoAPI.expect-Tuple{AbstractBlock,%20DensityMatrix})
in the `Yao.jl` documentation for the details.

## Sampling
Samples in the computational basis can be obtained with the same functions as 
in `Yao.jl`:

```@docs
measure(::MajoranaReg)
measure!(::MajoranaReg)
```

## Fidelities and bitstring probabilities
Fidelities between two `MajoranaReg`s can be computed the same way as in 
`Yao.jl`:

```@docs
fidelity(::MajoranaReg, ::MajoranaReg)
```

Under the hood this function uses a non-exported function to compute the 
probability of obtaining a given `bit_string` as an output when measuring in 
the computational basis:

```@docs
FLOYao.bitstring_probability
```

## GPU support
To enable the `FLOYaoCUDAExt` package extension you need to add [`CUDA.jl`](https://github.com/JuliaGPU/CUDA.jl), [`KernelAbstractions.jl`](https://github.com/JuliaGPU/KernelAbstractions.jl)
and [`ExponentialUtilities.jl`](https://github.com/SciML/ExponentialUtilities.jl)
to your project via
```julia-repl
julia> ]
(pgk) > add CUDA, KernelAbstractions, ExponentialUtilities
```
The reason that adding `CUDA.jl` to the project doesn't suffice is that 
[package extensions can't have dependencies](https://github.com/JuliaLang/Pkg.jl/issues/3641).

Once that is done, GPU support implemented very similar to the way it works in
`Yao.jl`. Array creation on the GPU is facilitated with the following functions:
```@docs
FLOYao.cuzero_state
FLOYao.cuone_state
FLOYao.cuproduct_state
```
which like [`FLOYao.one_state`](@ref) are not exported to not clash with 
the corresponding functions in `Yao.jl`. 

Upload and download of `MajoranaReg`s to and from the GPU is facilitated by
the `cu` and `cpu` functions
```@docs
FLOYao.cu
FLOYao.cpu
```
like in `Yao.jl`.

A quick example from the tests is shown below:
```julia
nq = 4
circuit = chain(nq)

θ = π/8
xxg = kron(nq, 1 => X, 2 => Y)
rg = rot(xxg, θ)
push!(circuit, rg)
push!(circuit, put(nq, 3=>Rz(0.5)))
push!(circuit, rg)

θ = π/5
xxg = kron(nq, 2 => X, 3 => Z, 4 => Y)
rg = rot(xxg, θ)
rz = put(nq, 3 => Rz(θ))
push!(circuit, rg)
push!(circuit, rz)

ham = put(nq, 1=>Z) + 2kron(nq, 1=>X, 2=>Z, 3=>Z, 4=>X) + 3.05f0put(nq, 2=>Z)
gpureg = FLOYao.cuzero_state(nq)
cpureg = FLOYao.zero_state(nq)
gpureg |> put(nq, 2=>X)
cpureg |> put(nq, 2=>X)


gpueval = expect(ham, gpureg |> circuit)
cpueval = expect(ham, cpureg |> circuit)
@assert meval ≈ aeval
@assert cpureg ≈ gpureg |> cpu
```

## Non-exported functions
The following functions are not exported and not really needed for the 
functionality of `FLOYao`, but can be useful for debugging.

```@docs
FLOYao.paulibasis2qubitop
FLOYao.majoranasquares2qubitbasis
FLOYao.qubit2paulibasis
FLOYao.qubit2majoranaevolution
FLOYao.paulibasis2majoranasquares
FLOYao.yaoham2majoranasquares
FLOYao.majorana2arrayreg
FLOYao.random_orthogonal_matrix
```
