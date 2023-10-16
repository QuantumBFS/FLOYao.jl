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
