# FLOYao.jl

A backend to efficiently simulate fermionic linear optics circuits (FLO) [Yao.jl](https://github.com/QuantumBFS/Yao.jlhttps://github.com/QuantumBFS/Yao.jl) based on [Classical simulation of noninteracting-fermion quantum circuits](https://arxiv.org/abs/quant-ph/0108010).

## Contents
 - [Installation](#Installation)
 - [Basic usage](#Basic-usage)
 - [Example: Transverse field Ising model](#Example:-Transverse-Field-Ising-model)
 - [Adding support for your own gates](#Adding-support-for-your-own-gates)
 - [Background: Fermionic linear optics circuits](#Background:-Fermionic-linear-optics-circuits)
 - [Known restrictions](#Known-restrictions)

## Installation
To install `FLOYao` open the julia REPL and type the following:
```julia
pkg> add FLOYao
```

## Basic usage
The heart of `FLOYao` is the `MajoranaReg` register type, which efficiently represents a state initialized as 
```julia
ψ = zero_state(nq)
ψ |> kron(nq, q => X for q in 1:k) 
```
and with a FLO circuit applied to it. Or expressed in the fermionic language:
$$
    \psi = \prod_{i=1}^k c_i^\dagger |\Omega\rangle 
$$

First import `Yao` and `FLOYao`:


```julia
using Yao, FLOYao
```

then build a (here somewhat arbitrary) circuit consisting only of [FLO gates](#Background:-Fermionic-linear-optics-circuits)


```julia
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

circuit
```




    nqubits: 4
    chain
    ├─ put on (3)
    │  └─ rot(Z, 0.5)
    ├─ rot(nqubits: 4
    kron
    ├─ 1=>X
    └─ 2=>X, 0.39269908169872414)
    ├─ rot(nqubits: 4
    kron
    ├─ 2=>X
    ├─ 3=>Z
    └─ 4=>X, 0.39269908169872414)
    ├─ put on (3)
    │  └─ rot(Z, 0.5)
    └─ put on (1)
       └─ Z




and define an observable that is a sum of squares of Majorana operators


```julia
hamiltonian = xxg1 + xxg2 + xxg3 + kron(nq, 2=>Z) + kron(nq, 3=>Z)
```




    nqubits: 4
    +
    ├─ kron
    │  ├─ 1=>X
    │  └─ 2=>X
    ├─ kron
    │  ├─ 2=>X
    │  ├─ 3=>Z
    │  └─ 4=>X
    ├─ kron
    │  ├─ 2=>X
    │  └─ 3=>X
    ├─ kron
    │  └─ 2=>Z
    └─ kron
       └─ 3=>Z




and finally create a register in the computational zero state via


```julia
mreg = MajoranaReg(nq, 0) # this is equivalent to `zero_state)`
```




    MajoranaReg{Float64} with 4 qubits:
    8×8 Matrix{Float64}:
     1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
     0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0
     0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0
     0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0
     0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0
     0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0



Applying the circuit to the register works then exactly the same way as for a normal Array register:


```julia
apply(mreg, circuit)
```




    MajoranaReg{Float64} with 4 qubits:
    8×8 Matrix{Float64}:
     -1.0  -0.0       -0.0       -0.0       -0.0       -0.0       -0.0       -0.0
     -0.0  -0.92388    0.382683  -0.0       -0.0       -0.0       -0.0       -0.0
      0.0   0.382683   0.92388    0.0        0.0        0.0        0.0        0.0
      0.0   0.0        0.0        0.92388    0.0        0.0        0.382683   0.0
      0.0   0.0        0.0        0.0        0.540302   0.841471   0.0        0.0
      0.0   0.0        0.0        0.0       -0.841471   0.540302   0.0        0.0
      0.0   0.0        0.0       -0.382683   0.0        0.0        0.92388    0.0
      0.0   0.0        0.0        0.0        0.0        0.0        0.0        1.0



and the same goes for expectation values of observables


```julia
expval = expect(hamiltonian, mreg => circuit)
```




    1.8535533905932737



or even gradients of these expectation values with respect to the circuit parameters


```julia
inδ, paramsδ = expect'(hamiltonian, mreg => circuit)
```




    MajoranaReg{Float64}(4, 0) => [0.0, -0.3535533905932738, -0.3535533905932738, 0.0]



and just to check that this is all consistent with running a full wavefunction simulation, can apply the same circuit on an `ArrayReg` and compare the expectation values and gradients


```julia
areg = zero_state(nq)
expval_full = expect(hamiltonian, areg => circuit)
expval_full ≈ expval
```




    true




```julia
inδ_full, paramsδ_full = expect'(hamiltonian, areg => circuit)
paramsδ ≈ paramsδ_full
```




    true



## Example: Transverse Field Ising model

For a more realistic use case, we have a look at VQE for the Transverse Field Ising model on a line whose Hamiltonian is given as 
$$
    H = J ∑_i^{L-1} X_i X_{i+1} + h ∑_i^L Z_i = U + T.
$$
and as Ansatz circuits we use the Hamiltonian Variational Ansatz
$$
    U(\vec θ) = ∏_i^p e^{-iθ_{i,U} U} e^{-iθ_{i,T} T} 
$$
with the initial state being the groundstate of the TFIM at $J = 0$, $|ψ_i⟩ = |0⋯0⟩$


```julia
# note that this is far beyond what is possible with a full wavefunction simulation
L = 100 
J = 0.5
h = 1.
p = 10

U = map(1:L-1) do i
    J * kron(L, i => X, i+1 => X)
end |> sum

T = map(1:L) do i
    h * kron(L, i => Z)
end |> sum

hamiltonian = T + U

circuit = chain(L)
for _ in 1:p
    for i in 1:L-1
        push!(circuit, rot(kron(L, i => X, i+1 => X), 0.))
    end
    for i in 1:L
        push!(circuit, put(L, i => Rz(0.)))
    end
end
nparams = nparameters(circuit)
dispatch!(circuit, rand(nparams) ./ 100)

ψ_i = MajoranaReg(L, 0);
```

now that we defined the hamiltonian, the ansatz circuit and the initial state we can perform
simple gradient descent on the energy expectation value to find an approximation to the
groundstate of $H$:


```julia
iterations = 100
γ = 2e-2

for i in 1:iterations
    _, grad = expect'(hamiltonian, ψ_i => circuit)
    dispatch!(-, circuit, γ * grad)
    println("Iteration $i, energy = $(expect(hamiltonian, ψ_i => circuit))")
end
```

    Iteration 1, energy = 99.72350957802038
    Iteration 2, energy = 99.40506967996363
    Iteration 3, energy = 98.76202821131267
    Iteration 4, energy = 97.46432218806858
    Iteration 5, energy = 94.8760933406418
    Iteration 6, energy = 89.85744274997147
    Iteration 7, energy = 80.66495993897014
    Iteration 8, energy = 65.52613565874718
    Iteration 9, energy = 44.706277096378145
    Iteration 10, energy = 22.702704698532308
    Iteration 11, energy = 5.585683940284514
    Iteration 12, energy = -4.957053736053638
    Iteration 13, energy = -11.51507490034233
    Iteration 14, energy = -16.938644625891662
    Iteration 15, energy = -22.799858577757277
    Iteration 16, energy = -29.61683483203969
    Iteration 17, energy = -37.178504470985246
    Iteration 18, energy = -44.84203010830907
    Iteration 19, energy = -51.911554043113405
    Iteration 20, energy = -57.94965135648679
    Iteration 21, energy = -62.88351660910522
    Iteration 22, energy = -66.88000315608417
    Iteration 23, energy = -70.13254847926378
    Iteration 24, energy = -72.74980240802692
    Iteration 25, energy = -74.78032473977396
    Iteration 26, energy = -76.27475030334335
    Iteration 27, energy = -77.32330501398673
    Iteration 28, energy = -78.05730844801298
    Iteration 29, energy = -78.6061254120199
    Iteration 30, energy = -79.06288335380546
    Iteration 31, energy = -79.48723263720359
    Iteration 32, energy = -79.92171373509012
    Iteration 33, energy = -80.40414299235896
    Iteration 34, energy = -80.97238940544219
    Iteration 35, energy = -81.65966007011164
    Iteration 36, energy = -82.47664457868399
    Iteration 37, energy = -83.38312975005772
    Iteration 38, energy = -84.27426661772805
    Iteration 39, energy = -85.02115376679511
    Iteration 40, energy = -85.55215722570551
    Iteration 41, energy = -85.88655795558829
    Iteration 42, energy = -86.0896267140955
    Iteration 43, energy = -86.22010996048115
    Iteration 44, energy = -86.31394312054839
    Iteration 45, energy = -86.38945436261524
    Iteration 46, energy = -86.45528566045766
    Iteration 47, energy = -86.51545669675187
    Iteration 48, energy = -86.57190690082466
    Iteration 49, energy = -86.6256492624962
    Iteration 50, energy = -86.67727270715429
    Iteration 51, energy = -86.72715905463484
    Iteration 52, energy = -86.77557842615364
    Iteration 53, energy = -86.82273313316587
    Iteration 54, energy = -86.86877942827263
    Iteration 55, energy = -86.91383943206021
    Iteration 56, energy = -86.9580084699872
    Iteration 57, energy = -87.00136011072875
    Iteration 58, energy = -87.04394995570384
    Iteration 59, energy = -87.08581868950496
    Iteration 60, energy = -87.12699465449951
    Iteration 61, energy = -87.1674960927351
    Iteration 62, energy = -87.20733313511003
    Iteration 63, energy = -87.2465095821556
    Iteration 64, energy = -87.28502449976742
    Iteration 65, energy = -87.32287364077607
    Iteration 66, energy = -87.36005069629678
    Iteration 67, energy = -87.39654837760624
    Iteration 68, energy = -87.43235932877198
    Iteration 69, energy = -87.46747687161835
    Iteration 70, energy = -87.50189558723514
    Iteration 71, energy = -87.53561174157423
    Iteration 72, energy = -87.56862356625369
    Iteration 73, energy = -87.60093140904806
    Iteration 74, energy = -87.63253777131091
    Iteration 75, energy = -87.66344725147033
    Iteration 76, energy = -87.69366641458836
    Iteration 77, energy = -87.7232036077549
    Iteration 78, energy = -87.75206873988294
    Iteration 79, energy = -87.7802730424737
    Iteration 80, energy = -87.80782882538476
    Iteration 81, energy = -87.83474923882456
    Iteration 82, energy = -87.8610480499644
    Iteration 83, energy = -87.88673943991061
    Iteration 84, energy = -87.91183782444902
    Iteration 85, energy = -87.93635770005173
    Iteration 86, energy = -87.96031351513375
    Iteration 87, energy = -87.98371956545809
    Iteration 88, energy = -88.00658991185603
    Iteration 89, energy = -88.02893831799855
    Iteration 90, energy = -88.05077820575987
    Iteration 91, energy = -88.07212262568603
    Iteration 92, energy = -88.09298424017994
    Iteration 93, energy = -88.11337531717759
    Iteration 94, energy = -88.13330773230405
    Iteration 95, energy = -88.15279297771876
    Iteration 96, energy = -88.17184217608879
    Iteration 97, energy = -88.19046609834005
    Iteration 98, energy = -88.2086751840374
    Iteration 99, energy = -88.22647956342138
    Iteration 100, energy = -88.24388908028976


## Adding support for your own gates

Natively, the only FLO gates that come already shipped with `Yao.jl` are gates of
the form `kron(nq, i => σ1, i+1 => Z, ⋯, i-1 => Z, j => σ2)` and 
`rot(kron(nq, i => σ1, i+1 => Z, ⋯, i-1 => Z, j => σ2), θ)` where `σ1, σ2 ∈ [X,Y]`. But there are many more FLO gates, one being for example the $FSWAP$ gates which swaps to qubits while making sure to preserve the fermionic commutation relations


```julia
@const_gate FSWAP::ComplexF64 = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 -1]
```

If a gate defines a matrix representation, as we just did for the `FSWAP`gate, `FLOYao` supports them out of the box by manually checking if they are a FLO gate and then computing its matrix representation in the Majorana basis. But this method is fairly slow, compared to directly implementing `unsafe_apply!(::MajoranaReg, ::YourBlock)` and  `instruct!(::MajoranaReg, ::YourBlock)` and will warn you accordingly


```julia
nq = 4
fswap = put(nq, (1, 2) => FSWAP)
mreg = MajoranaReg(nq, 1)
mreg |> fswap
```

    ┌ Warning: Calling manual instruct!(MajoranaReg{Float64}(4, 1), ComplexF64[1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im -1.0 + 0.0im], (1, 2)).
    │ You can greatly speed up your FLO gates by exactly implementing unsafe_apply!()
    │ and instruct!() for them. See FLOYao/src/instruct.jl for how to do that.
    └ @ FLOYao /home/yc20910/PhD/Work/code/FLOYao/src/instruct.jl:49





    MajoranaReg{Float64} with 4 qubits:
    8×8 Matrix{Float64}:
     -2.40901e-16  -2.94663e-16  -1.0           2.35127e-16  0.0  0.0  0.0  0.0
      3.57633e-16  -4.32975e-16  -1.0529e-16   -1.0          0.0  0.0  0.0  0.0
     -1.0           1.20422e-17   2.61492e-16  -2.81188e-16  0.0  0.0  0.0  0.0
     -7.05281e-17  -1.0           2.65282e-16   3.00191e-16  0.0  0.0  0.0  0.0
      0.0           0.0           0.0           0.0          1.0  0.0  0.0  0.0
      0.0           0.0           0.0           0.0          0.0  1.0  0.0  0.0
      0.0           0.0           0.0           0.0          0.0  0.0  1.0  0.0
      0.0           0.0           0.0           0.0          0.0  0.0  0.0  1.0



now before we fix these warnings, let's see how long the current implementation takes:


```julia
using BenchmarkTools
using Suppressor # we don't want to get all the warnings when benchmarking
@benchmark @suppress apply!($mreg, $fswap)
```




    BenchmarkTools.Trial: 6960 samples with 1 evaluation.
     Range (min … max):  639.587 μs …   4.007 ms  ┊ GC (min … max): 0.00% … 70.64%
     Time  (median):     677.467 μs               ┊ GC (median):    0.00%
     Time  (mean ± σ):   713.191 μs ± 217.746 μs  ┊ GC (mean ± σ):  2.39% ±  6.42%
    
      ▄▅▆▇█▇▅▅▄▃▂▂▂▁▁▁▁                                             ▂
      █████████████████▆▇▆▆▆▆▅▅▅▅▆▄▅▄▅▅▄▄▄▇▆▆▆▇▆▅▆▄▅▄▅▅▅▅▅▅▄▄▆▅▄▆▃▄ █
      640 μs        Histogram: log(frequency) by time       1.15 ms <
    
     Memory estimate: 339.34 KiB, allocs estimate: 501.



To find out what the matrix representation of the `FSWAP` gate in the Majorana basis is, it is easiest to retrace what is happening inside `instruct!(::MajoranaReg, ::AbstractMatrix, locs)`


```julia
W = FLOYao.qubit2majoranaevolution(Matrix(fswap.content), fswap.locs)
```




    4×4 Matrix{Float64}:
     -2.40901e-16  -2.94663e-16  -1.0           2.35127e-16
      3.57633e-16  -4.32975e-16  -1.0529e-16   -1.0
     -1.0           1.20422e-17   2.61492e-16  -2.81188e-16
     -7.05281e-17  -1.0           2.65282e-16   3.00191e-16




```julia
matlocs = 2*(fswap.locs[1]-1)+1:2(fswap.locs[end])
```




    1:4



this matrix gets left-multiplied to the columns `1:4` in the last line of `FLOYao.majorana_unsafe_apply!(::MajoranaReg, ::PutBlock)`. So we can instead implement the action of an `FSWAP` gate on a `MajoranaReg` directly as follows:


```julia
function YaoBlocks.unsafe_apply!(reg::MajoranaReg, b::PutBlock{2,2,FSWAPGate})
    FLOYao.areconsecutive(b.locs) || throw(NonFLOException("FSWAP must act on consecutive qubits"))
    instruct!(reg, Val(:FSWAP), b.locs)
end

function Yao.instruct!(reg::MajoranaReg, ::Val{:FSWAP}, locs::Tuple)
    i1, i2 = locs
    ψ1, ψ2 = reg.state[2i1-1,:], reg.state[2i1,:]
    ψ3, ψ4 = reg.state[2i2-1,:], reg.state[2i2,:]
    reg.state[2i1-1,:] .=  .-ψ3
    reg.state[2i1,:] .=  .-ψ4
    reg.state[2i2-1,:] .=  .-ψ1
    reg.state[2i2,:] .=  .-ψ2
    return reg
end
```


```julia
@benchmark apply!($mreg, $fswap)
```




    BenchmarkTools.Trial: 10000 samples with 709 evaluations.
     Range (min … max):  175.340 ns …   3.170 μs  ┊ GC (min … max): 0.00% … 94.01%
     Time  (median):     186.623 ns               ┊ GC (median):    0.00%
     Time  (mean ± σ):   208.873 ns ± 179.849 ns  ┊ GC (mean ± σ):  7.50% ±  8.05%
    
      ▃█▅▆▇▅▅▄▃▂▁▁▁                                                 ▂
      ██████████████▇▆▆▄▆▅▅▅▄▄▅▄▃▄▁▄▁▁▁▃▁▃▃▁▃▄▁▁▁▁▁▅▇▇█▆█████▇█▇▆▅▅ █
      175 ns        Histogram: log(frequency) by time        370 ns <
    
     Memory estimate: 512 bytes, allocs estimate: 4.



Now, that is quite a significant speedup!

## Background: Fermionic linear optics circuits

This section is more here, to fix the convention of
[Jordan-Wigner transform](https://en.wikipedia.org/wiki/Jordan%E2%80%93Wigner_transformation)
and [Majorana operators](https://en.wikipedia.org/wiki/Majorana_fermion) that we use here,
and less to explain the full theory behind those. For the latter, we, once again, recommend [Classical simulation of noninteracting-fermion quantum circuits](https://arxiv.org/abs/quant-ph/0108010).

We define the Majorana operators $γ_i$ via 
$$
    γ_{2i-1} = ∏_{j=1}^{i-1} (-Z_j) X_i
$$
and
$$
    γ_{2i} = ∏_{j=1}^{i-1} (-Z_j) Y_i.
$$
This implies the normal fermionic creation and annihilation operators are given by
$$
    c_j = \frac{1}{2} (γ_{2j-1} + iγ_{2j})
    \quad \textrm{and} \quad
    c_j^† = \frac{1}{2} (γ_{2j-1} - iγ_{2j})
$$
and products of two Majorana operators are of the form
$$
    σ_i \left(∏_{i<j<k} Z_k \right) σ_k
    \quad \textrm{or} \quad
    Z_i
$$
with $σ_i, σ_k ∈ \{X, Y\}$.

If a unitary is of the form 
$$
    U = e^{-iθH}
$$
with 
$$
    H = \frac{i}{4} \sum_{i,j} H^{ij} γ_i γ_j
$$
it turns out that it takes Majorana operators to linear combinations of 
Majorana operators under conjugation, i.e.
$$
    U γ_i U^† = R_i^j γ_j
$$
with some $R ∈ SO(2n)$. The action of FLO circuits on input states of the form 
$$
    |ψ_i⟩ = c_k^† ⋯ c_1^† |0 ⋯ 0⟩ =  γ_{2k-1} ⋯ γ_{2k-1} |0 ⋯ 0⟩
$$
can thus be efficiently classically simluated.

Furthermore, it is also fairly straightforward to compute expectation values since we can 
evolve the Hamiltonian first in the Heisenber picture to
$$
    UHU^† = \frac{i}{4} R^{m}_{i} R^{n}_{j} H^{ij} γ_{m} γ_{n} 
           =: \frac{i}{4} \tilde H^{mn} γ_{m} γ_{n}.
$$

For the expectation value this makes then
$$
\begin{aligned}
    ⟨ψ|UHU^†|ψ⟩ &= \frac{i}{4} ∑_i \tilde H^{ii} ⟨Ω|γ_{1} ⋯ γ_{2k-1} γ_{i} γ_{i} γ_{2k-1} ⋯  γ_{1}|Ω⟩ \\
                &= \frac{i}{2} ∑_{i ≤ k} \tilde H^{2i-1,2i} ⟨Ω| γ_{2i-1} γ_{2i-1} γ_{2i} γ_{2i-1}|Ω⟩
                   + \frac{i}{2} ∑_{i > k} \tilde H^{2i-1,2i} ⟨Ω|γ_{2i-1} γ_{2i}|Ω⟩ \\
                &= \frac{1}{2} ∑_{i ≤ k} \tilde H^{2i-1,2i} - \frac{1}{2} ∑_{i>k} \tilde H^{2i-1,2i} \\
                &= \frac{1}{2} ∑_{i ≤ k} R^{2i-1}_{m} R^{2i}_{n} H^{mn} 
                   - \frac{1}{2} ∑_{i>k} R^{2i-1}_{m} R^{2i}_{n} H^{mn} \\
\end{aligned}
$$
where we used that the fermionic vacuum $|Ω⟩$ is the zero state.


## Known restrictions

###  Expectation values of higher order observables
So far, `FLOYao` only supports expectation values of observables that are sums of squares of 
Majorana operators. But in general, one can use [Wick's theorem](https://en.wikipedia.org/wiki/Wick%27s_theorem) to calculate the expectation values of expressions of the form 
$$
    ⟨ψ_i|O|ψ_i⟩
$$
where $|ψ_i⟩$ is a computational basis state and $O$ a string of Majorana operators and thus using linearity of the expectation value it is possible to efficiently calculate the expectation value of any observable that can be expanded in a sum of polynomially (in the number of qubits) many products of Majorana operators. (See also [Classical simulation of noninteracting-fermion quantum circuits](https://arxiv.org/abs/quant-ph/0108010) again for details). If you need expectation values of higher order (in the number of Majorana operators involved) observables, feel free to open a PR!

### "Hidden" FLO circuits
`Yao.jl` allows to express the same circuit in different forms. E.g. 
```julia
    chain(nq, put(1 => X), put(2 => X))
```
and
```
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

If you run into a case that is a FLO circuit / gate but not recognised as such please open an issue or even pull request.

