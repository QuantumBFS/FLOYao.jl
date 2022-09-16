# Known restrictions

##  Expectation values of higher order observables
So far, `FLOYao` only supports expectation values of observables that are sums of squares of 
Majorana operators. But in general, one can use [Wick's theorem](https://en.wikipedia.org/wiki/Wick%27s_theorem)
to calculate the expectation values of expressions of the form 
```math
    ⟨O⟩ = ⟨\Omega |U^† O U|\Omega ⟩
```
where $|Ω⟩ = |0 ⋯ 0⟩$ is the all zero stat, $O$ a string of Majorana operators
and $U$ a FLO unitary. Thus, using linearity of the expectation value, it is
possible to efficiently calculate the expectation value of any observable that
can be expanded in a sum of polynomially (in the number of qubits) many
products of Majorana operators. (See also [Classical simulation of noninteracting-fermion quantum circuits](https://arxiv.org/abs/quant-ph/0108010) again for details).
If you need expectation values of higher order (in the number of Majorana
operators involved) observables, feel free to open a pull request!

## "Hidden" FLO circuits
`Yao.jl` allows to express the same circuit in different forms. E.g. 
```julia
chain(nq, put(1 => X), put(2 => X))
```
and
```julia
kron(nq, 1 => X, 2 => X)
```
represent the same FLO circuit, but only the latter will be recognised as such. Similarly 
```julia
kron(nq, 1 => X, 2 => X, 3 => X, 4 => X)
```
and
```julia
chain(nq, kron(1 => X, 2 => X), kron(3 => X, 4 => X))
```
represent the same FLO circuit but only the latter will be recognised as such. This is because

 - We don't check if a whole `ChainBlock` is a FLO circuit, even if its single
   gates are not. Instead a `ChainBlock` is applied gate by gate, each of which
   has to be a FLO gate.
 - For `KronBlock`s we first try if each factor is a FLO gate and then if the whole
   block altogether is of the form  `kron(nq, i => σ1, i+1 => Z, ⋯, i-1 => Z, j => σ2)`
   with `σ1, σ2 ∈ [X, Y]`.

If you run into a case that is a FLO circuit / gate but not recognised as such
please open an issue or even pull request.

