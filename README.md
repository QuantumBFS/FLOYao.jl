# FLOYao.jl

[![CI][ci-img]][ci-url]
[![codecov][codecov-img]](codecov-url)
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

A brief introduction to fermionic linear optics circuits is found in the 
[Documentation](docs-stable-url) and a more in-depth introduction in e.g. the two papers linked above.


[ci-img]: https://github.com/QuantumBFS/FLOYao.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/QuantumBFS/FLOYao.jl/actions
[codecov-img]: https://codecov.io/gh/QuantumBFS/FLOYao.jl/branch/master/graph/badge.svg?token=U604BQGRV1
[codecov-url]: https://codecov.io/gh/QuantumBFS/FLOYao.jl
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://QuantumBFS.github.io/FLOYao.jl/dev/
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://QuantumBFS.github.io/FLOYao.jl/stable
