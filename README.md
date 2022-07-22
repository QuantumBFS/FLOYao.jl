# FLOYao.jl

A backend to efficiently simulate fermionic linear optics (FLO) circuits in [Yao.jl](https://github.com/QuantumBFS/Yao.jlhttps://github.com/QuantumBFS/Yao.jl) based on [Classical simulation of noninteracting-fermion quantum circuits](https://arxiv.org/abs/quant-ph/0108010).

**Note**    
The markdown version of this README is automatically generated from `README.ipynb` and some 
of the hpyerlinks and math doesn't seem to play that well with githubs markdown parser. So you might want to have a look at the `README.ipynb` with jupyter instead for a better reading experience.

## Contents
 - [Installation](#Installation)
 - [Basic usage](#Basic-usage)
 - [List of supported gates](#List-of-supported-gates)
 - [Example: Transverse field Ising model](#Example:-Transverse-Field-Ising-model)
 - [Adding support for your own gates](#Adding-support-for-your-own-gates)
 - [Background: Fermionic linear optics circuits](#Background:-Fermionic-linear-optics-circuits)
 - [Known restrictions](#Known-restrictions)

## Installation

 > By the time this will be in the julia general registry, it should be as easy as opening
 > a julia REPL and typing
 >  ```julia
 >   pkg> add FLOYao
 >  ```

But for now I suggest the following way: First open a standard terminal and type 
```bash
cd your_favorite_folder_for_code
git clone git@github.com:PhaseCraft/FLOYao.jl
```
and then open a julia REPL and type
```
pkg> dev your_favorite_folder_for_code/FLOYao.jl
```
which should make `FLOYao.jl` discoverable for your standard julia installation. 
Under linux the standard folder for julia packages under development is `/home/username/.julia/dev`

## Basic usage
The heart of `FLOYao` is the `MajoranaReg` register type, which efficiently represents a state that is a [FLO unitary](#Background:-Fermionic-linear-optics-circuits) applied to the vacuum state $|0⋯0⟩$

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
mreg = MajoranaReg(nq)
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



Applying the circuit to the register works then exactly the same way as for a normal `ArrayReg` register:


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




    MajoranaReg{Float64}(4) => [0.0, -0.3535533905932738, -0.3535533905932738, 0.0]



Just to check that this is all consistent with running a full wavefunction simulation, we can apply the same circuit on an `ArrayReg` and compare the expectation values and gradients


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



## List of supported gates

The following gates are FLO gates and supported by `FLOYao.jl`:

|  Gate | Comment |
|-------|---------|
|  `XGate`   | Together with `Y` the only gate that does not preserve fermionic parity |
|  `YGate`   |   See above  |
|  `ZGate`   |         |
|  `RotationGate{⋯,⋯,YGate}`  | The only single qubit rotation gate since $R_x(θ)γ_i R_x(-θ)$ is not a linear combination of Majorana operators for all Majorana operators. Similar for $R_y$ |
| `PauliKronBlock` | A kronecker product of Pauli operators s.t. that first and last operator are either $X$ or $Y$ and all operators in between are $Z$.  |
| `RotationGate{⋯,⋯,PauliKronBlock}` | A rotation gate with generator as in the last line. |
| `AbstractMatrix` | Unless the gate type is already explicitely implemented or know to not be a FLO gate, `FLOYao` will try to automatically convert the gate matrix in the qubit basis to a matrix in the Majorana basis. But note that this is fairly slow (although still poly-time instead of exp-time) |


If you want to add support to your own gates, read [this section](#Adding-support-for-your-own-gates) to learn how to do that.

## Example: VQE for the Transverse Field Ising model

For a more realistic use case, we have a look at VQE for the Transverse Field Ising model on a line whose Hamiltonian is given as 
$$
    H = J ∑_i^{L-1} X_i X_{i+1} + h ∑_i^L Z_i = U + T.
$$
As Ansatz circuits we use the Hamiltonian Variational Ansatz
$$
    U(\vec θ) = ∏_i^p e^{-iθ_{i,U} U} e^{-iθ_{i,T} T} 
$$
with the initial state being the groundstate of the TFIM at $J = 0$, so $|ψ_i⟩ = |0⋯0⟩$


```julia
L = 100 # this is far beyond what is possible with a full wavefunction simulation
J = 0.5
h = -1.
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

ψ_i = MajoranaReg(L);
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

    Iteration 1, energy = -99.88490767843159
    Iteration 2, energy = -99.97441634919842
    Iteration 3, energy = -100.0056763860236
    Iteration 4, energy = -100.01859116906469
    Iteration 5, energy = -100.02599635617956
    Iteration 6, energy = -100.03205626486043
    Iteration 7, energy = -100.03817176257458
    Iteration 8, energy = -100.04487079247842
    Iteration 9, energy = -100.05241364153315
    Iteration 10, energy = -100.0609878071107
    Iteration 11, energy = -100.07077199153714
    Iteration 12, energy = -100.08195848715462
    Iteration 13, energy = -100.09476245263055
    Iteration 14, energy = -100.10942717398012
    Iteration 15, energy = -100.12622819597735
    Iteration 16, energy = -100.1454772116763
    Iteration 17, energy = -100.16752593869434
    Iteration 18, energy = -100.19276997288964
    Iteration 19, energy = -100.2216524987739
    Iteration 20, energy = -100.25466765628717
    Iteration 21, energy = -100.29236328419954
    Iteration 22, energy = -100.33534267125965
    Iteration 23, energy = -100.38426484506415
    Iteration 24, energy = -100.43984281764033
    Iteration 25, energy = -100.50283909225735
    Iteration 26, energy = -100.57405762939587
    Iteration 27, energy = -100.65433138893495
    Iteration 28, energy = -100.74450453582841
    Iteration 29, energy = -100.84540845122724
    Iteration 30, energy = -100.95783087041501
    Iteration 31, energy = -101.08247781644765
    Iteration 32, energy = -101.21992855313191
    Iteration 33, energy = -101.37058456596847
    Iteration 34, energy = -101.53461458592093
    Iteration 35, energy = -101.71189883936002
    Iteration 36, energy = -101.901976913864
    Iteration 37, energy = -102.10400467832382
    Iteration 38, energy = -102.31672633328924
    Iteration 39, energy = -102.53846761974292
    Iteration 40, energy = -102.76715525120282
    Iteration 41, energy = -103.00036564869174
    Iteration 42, energy = -103.2354031419488
    Iteration 43, energy = -103.4694042807764
    Iteration 44, energy = -103.69946131713009
    Iteration 45, energy = -103.92275492642514
    Iteration 46, energy = -104.13668445415813
    Iteration 47, energy = -104.3389838195596
    Iteration 48, energy = -104.52781277778378
    Iteration 49, energy = -104.70181627432959
    Iteration 50, energy = -104.8601485638006
    Iteration 51, energy = -105.00246290042799
    Iteration 52, energy = -105.12887124926024
    Iteration 53, energy = -105.23988109223869
    Iteration 54, energy = -105.33631774759996
    Iteration 55, energy = -105.41924069062468
    Iteration 56, energy = -105.48986138550025
    Iteration 57, energy = -105.5494684692163
    Iteration 58, energy = -105.59936415865161
    Iteration 59, energy = -105.64081382274011
    Iteration 60, energy = -105.67500901690448
    Iteration 61, energy = -105.70304304828547
    Iteration 62, energy = -105.72589735822787
    Iteration 63, energy = -105.74443662995212
    Iteration 64, energy = -105.75941046942594
    Iteration 65, energy = -105.77145966743821
    Iteration 66, energy = -105.78112533829226
    Iteration 67, energy = -105.7888595704421
    Iteration 68, energy = -105.79503656333311
    Iteration 69, energy = -105.79996352961251
    Iteration 70, energy = -105.8038908962753
    Iteration 71, energy = -105.80702153772377
    Iteration 72, energy = -105.80951892120831
    Iteration 73, energy = -105.8115141478428
    Iteration 74, energy = -105.81311193899111
    Iteration 75, energy = -105.8143956569036
    Iteration 76, energy = -105.81543146764064
    Iteration 77, energy = -105.81627175984562
    Iteration 78, energy = -105.81695792968344
    Iteration 79, energy = -105.81752263387214
    Iteration 80, energy = -105.81799160173874
    Iteration 81, energy = -105.81838508535368
    Iteration 82, energy = -105.818719015115
    Iteration 83, energy = -105.81900591730428
    Iteration 84, energy = -105.81925564043273
    Iteration 85, energy = -105.81947592876408
    Iteration 86, energy = -105.81967287421094
    Iteration 87, energy = -105.81985127178436
    Iteration 88, energy = -105.82001489879488
    Iteration 89, energy = -105.82016673392603
    Iteration 90, energy = -105.82030912899383
    Iteration 91, energy = -105.82044394353889
    Iteration 92, energy = -105.82057265026401
    Iteration 93, energy = -105.82069641762584
    Iteration 94, energy = -105.82081617454064
    Iteration 95, energy = -105.82093266109534
    Iteration 96, energy = -105.82104646831016
    Iteration 97, energy = -105.82115806934105
    Iteration 98, energy = -105.82126784398505
    Iteration 99, energy = -105.82137609794717
    Iteration 100, energy = -105.82148307800819


## Adding support for your own gates

Natively, the only FLO gates that come already shipped with `Yao.jl` are the gates listed [here](#List-of-supported-gates). But there are many more FLO gates, one being for example the `FSWAP` gate which swaps to qubits while making sure to preserve the fermionic commutation relations


```julia
@const_gate FSWAP::ComplexF64 = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 -1]
```

If a gate defines a matrix representation, as we just did for the `FSWAP`gate, `FLOYao` supports them out of the box by manually checking if they are a FLO gate and then computing its matrix representation in the Majorana basis. But this method is fairly slow, compared to directly implementing `unsafe_apply!(::MajoranaReg, ::YourBlock)` and  `instruct!(::MajoranaReg, ::YourBlock)` and will warn you accordingly


```julia
nq = 4
fswap = put(nq, (1, 2) => FSWAP)
mreg = MajoranaReg(nq)
mreg |> put(nq, 2 => X)
mreg |> fswap
```

    ┌ Warning: Calling manual instruct!(MajoranaReg{Float64}(4), ComplexF64[1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im -1.0 + 0.0im], (1, 2)).
    │ You can greatly speed up your FLO gates by exactly implementing unsafe_apply!()
    │ and instruct!() for them. See FLOYao/src/instruct.jl for how to do that.
    └ @ FLOYao /home/yc20910/PhD/Work/code/FLOYao/src/instruct.jl:62





    MajoranaReg{Float64} with 4 qubits:
    8×8 Matrix{Float64}:
     -2.40901e-16  -2.94663e-16  -1.0          …   0.0   0.0   0.0   0.0
      3.57633e-16  -4.32975e-16  -1.0529e-16       0.0   0.0   0.0   0.0
     -1.0           1.20422e-17   2.61492e-16      0.0   0.0   0.0   0.0
     -7.05281e-17  -1.0           2.65282e-16      0.0   0.0   0.0   0.0
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




    BenchmarkTools.Trial: 7192 samples with 1 evaluation.
     Range (min … max):  646.570 μs …   3.542 ms  ┊ GC (min … max): 0.00% … 71.01%
     Time  (median):     658.311 μs               ┊ GC (median):    0.00%
     Time  (mean ± σ):   690.038 μs ± 194.914 μs  ┊ GC (mean ± σ):  2.22% ±  6.28%
    
      ▆█▆▄▃▃▂▁ ▁                                                    ▁
      ██████████▇▇▇▇▅▇▇▇▅▅▆▅▄▄▄▄▄▄▆▅▅▄▄▅▇▇█▇▇▇▆▅▅▃▄▃▃▃▃▁▁▁▁▁▁▃▁▅▆▇▆ █
      647 μs        Histogram: log(frequency) by time        1.1 ms <
    
     Memory estimate: 339.19 KiB, allocs estimate: 496.



To find out what the matrix representation of the `FSWAP` gate in the Majorana basis is, it is easiest to retrace what is happening inside `instruct!(::MajoranaReg, ::AbstractMatrix, locs)`
You can use


```julia
@which instruct!(mreg, mat(FSWAP), (1,2))
```




instruct!(reg::<b>MajoranaReg</b>, gate::<b>AbstractMatrix</b>, locs) in FLOYao at <a href="https://github.com/PhaseCraft/FLOYao.jl/tree/79ab50718cfa84dfdde4a4141af0d8beece24395//src/instruct.jl#L55" target="_blank">/home/yc20910/PhD/Work/code/FLOYao/src/instruct.jl:55</a>



to find the location of the corresponding code. Now let's copy-paste what we found there:


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
     Range (min … max):  170.410 ns …   2.536 μs  ┊ GC (min … max): 0.00% … 92.06%
     Time  (median):     181.900 ns               ┊ GC (median):    0.00%
     Time  (mean ± σ):   200.316 ns ± 151.778 ns  ┊ GC (mean ± σ):  6.67% ±  7.97%
    
       ▁▂█▆▅▆▆▄▃▂▁                                                  ▂
      ▇████████████▇▇▆▅▅▄▃▄▁▄▄▁▅▁▃▁▄▃▁▁▁▄▃▃▄▄▄▃▃▃▃▃▄▁▃▄▃▃▅▄▄▇▇▇▇▆▅▅ █
      170 ns        Histogram: log(frequency) by time        339 ns <
    
     Memory estimate: 512 bytes, allocs estimate: 4.



Now, that is quite a significant speedup!

## Background: Fermionic linear optics circuits

This section is here, more to fix the convention of
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

Any unitary that takes all Majorana operators to a linear combination of Majorana operators
under conjugation, i.e. that satisfies
$$
    U γ_i U^† = R_i^j γ_j
$$
with some $R ∈ O(2n)$ is a FLO unitary. In particular, if a unitary is of the form 
$$
    U = e^{-iθH}
$$
with 
$$
    H = \frac{i}{4} \sum_{i,j} H^{ij} γ_i γ_j
$$
it is a FLO unitary with even $R ∈ SO(2n)$.

But note, that not all FLO unitaries are of that form. For example $X_i$ is also a FLO 
gate since it either commutes or anti-commutes with all Majorana operators, but the associated
matrix $R$ always has determinant $-1$.

Calculating the expectation values of hamiltonians like the one above when evolving the 
vacuum state with FLO circuits is efficiently possible. First evolve the 
Hamiltonian in the Heisenber picture to
$$
    UHU^† = \frac{i}{4} R^{m}_{i} R^{n}_{j} H^{ij} γ_{m} γ_{n} 
           =: \frac{i}{4} \tilde H^{mn} γ_{m} γ_{n}.
$$
and then compute the expectation value
$$
\begin{aligned}
    ⟨ψ|UHU^†|ψ⟩ &= \frac{i}{4} \tilde H^{mn} ⟨Ω|γ_{m} γ_{n}|Ω⟩ \\
                &= - \frac{1}{2} ∑_{i} \tilde H^{2i-1,2i} \\
                &= - \frac{1}{2} ∑_{i>k} R^{2i-1}_{m} R^{2i}_{n} H^{mn} \\
\end{aligned}.
$$
From the first to second line one needs to carefully think which of the 
$⟨Ω|γ_{m} γ_{n}|Ω⟩$ are zero and which cancel each other out due to the anti-symmetry of $H^{mn}$.

## Known restrictions

###  Expectation values of higher order observables
So far, `FLOYao` only supports expectation values of observables that are sums of squares of 
Majorana operators. But in general, one can use [Wick's theorem](https://en.wikipedia.org/wiki/Wick%27s_theorem) to calculate the expectation values of expressions of the form 
$$
    ⟨ψ_i|O|ψ_i⟩
$$
where $|ψ_i⟩$ is a computational basis state and $O$ a string of Majorana operators and thus, using linearity of the expectation value, it is possible to efficiently calculate the expectation value of any observable that can be expanded in a sum of polynomially (in the number of qubits) many products of Majorana operators. (See also [Classical simulation of noninteracting-fermion quantum circuits](https://arxiv.org/abs/quant-ph/0108010) again for details). If you need expectation values of higher order (in the number of Majorana operators involved) observables, feel free to open a pull request!

### "Hidden" FLO circuits
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

If you run into a case that is a FLO circuit / gate but not recognised as such please open an issue or even pull request.

