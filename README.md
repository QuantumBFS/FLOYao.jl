# FLOYao.jl

[![CI][ci-img]][ci-url]
[![codecov][codecov-img]][codecov-url]
[![][docs-stable-img]][docs-stable-url]
[![][docs-dev-img]][docs-dev-url]


A [Yao.jl](https://github.com/QuantumBFS/Yao.jl) backend to efficiently simulated
fermionic linear optics (FLO) circuits in  based on
[Classical simulation of noninteracting-fermion quantum circuits](https://arxiv.org/abs/quant-ph/0108010)
and [Disorder-assisted error correction in Majorana chains](https://arxiv.org/abs/1108.3845).
FLO circuits are a class of quantum circuits that are closely related to
non-interacting fermions and can be efficiently simulated on classical
computers, similar to the way Clifford circuits can be efficiently classically
simulated, as is done in
[YaoClifford.jl](https://github.com/QuantumBFS/YaoClifford.jl).

The goal of `FLOYao` is that if you have code written in `Yao.jl` that only 
uses [FLO gates](https://quantumbfs.github.io/FLOYao.jl/stable/supported_gates/)
and other primitives that are efficiently simulatable in polynomial time and 
space, that you can simply replace your `AbstractArrayReg` with a `MajoranaReg`
and run exactly the same simulation, with the same code but exponentially faster.

A brief introduction to fermionic linear optics circuits is found in the 
[Documentation](docs-stable-url) and a more in-depth introduction in e.g. the two papers linked above.


## Installation
`FLOYao` can be simply installed from the REPL via

```jl-repl
pkg> add FLOYao
```

## Running circuits
First import `FLOYao` and `Yao`
```julia
using FLOYao, Yao
```
then build a (here somewhat arbitrary) circuit consisting only of [Supported gates](https://quantumbfs.github.io/FLOYao.jl/stable/background/)

```jldoctest quickstart; output=false
nq = 4
θ = π/8
circuit = chain(nq)

push!(circuit, put(nq, 3=>Rz(0.5)))

xxg1 = kron(nq, 1 => X, 2 => X)
rg = rot(xxg1, θ)
push!(circuit, rg)  

xxg2 = kron(nq, 2 => X, 3 => Z, 4 => X)
rg = rot(xxg2, θ)
push!(circuit, rg)  
push!(circuit, put(nq, 3=>Rz(0.5)))
push!(circuit, put(nq, 1=>Z))

xxg3 = kron(nq, 2 => X, 3 => X)
rg = rot(xxg3, θ)
push!(circuit, rg)
```

and create a FLO state, pipe it through the circuit and measure the result

```julia
FLOYao.zero_state(nq) |> circuit |> measure!
```

## Documentation
The documentation for the last release is [here][docs-stable-url] and the documentation
for the current development branch [here][docs-dev-url].

[ci-img]: https://github.com/QuantumBFS/FLOYao.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/QuantumBFS/FLOYao.jl/actions
[codecov-img]: https://codecov.io/gh/QuantumBFS/FLOYao.jl/branch/master/graph/badge.svg?token=U604BQGRV1
[codecov-url]: https://codecov.io/gh/QuantumBFS/FLOYao.jl
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://QuantumBFS.github.io/FLOYao.jl/dev/
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://QuantumBFS.github.io/FLOYao.jl/stable


