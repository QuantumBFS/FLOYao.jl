# FLOYao.jl

A backend to efficiently simulate fermionic linear optics (FLO) circuits in [Yao.jl](https://github.com/QuantumBFS/Yao.jlhttps://github.com/QuantumBFS/Yao.jl) based on [Classical simulation of noninteracting-fermion quantum circuits](https://arxiv.org/abs/quant-ph/0108010)
and [Disorder-assisted error correction in Majorana chains](https://arxiv.org/abs/1108.3845).

**Note**    
The markdown version of this README is automatically generated from `README.ipynb` and some 
of the hpyerlinks and math doesn't seem to play that well with githubs markdown parser. So you might want to have a look at `README.ipynb` with jupyter instead for a better reading experience.

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
mreg = FLOYao.zero_state(nq)
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



Other options to create registers are `FLOYao.product_state` and `FLOYao.zero_state_like`. `FLOYao.rand_state` is not implemented so far, because to do so I would need to rely on yet another dependency to create random orthogonal matrices.

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



It is also possible to sample from the state described by a `MajoranaReg`ister using `Yao`'s
`measure` or `measure!` functions. 


```julia
copyreg = copy(mreg) |> circuit
samples = measure(copyreg, nshots=10000)
samples[1:5]
```




    5-element Vector{DitStr{2, 4, BigInt}}:
     0000 ₍₂₎
     0000 ₍₂₎
     0000 ₍₂₎
     0000 ₍₂₎
     0000 ₍₂₎



and if you want to measure only a subset of qubits or change the ordering in
which the qubits are measured, this is also possible:


```julia
samples213 = measure(copyreg, [2, 1, 3],  nshots=10000)
samples213[1:5]
```




    5-element Vector{DitStr{2, 3, BigInt}}:
     000 ₍₂₎
     000 ₍₂₎
     000 ₍₂₎
     000 ₍₂₎
     000 ₍₂₎



Just to check that this is all consistent with running a full wavefunction simulation, we can apply the same circuit on an `ArrayReg` and compare the expectation values, gradients and samples


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




```julia
copyareg = copy(areg) |> circuit
samples_full = measure(copyareg, nshots=10000)
samples213_full = measure(copyareg, [2, 1, 3], nshots=10000);
```


```julia
using StatsBase
println("Comparing the full countmaps:")
println("")
display(countmap(samples))
display(countmap(samples_full))
println("-----------------------------")

println("")
println("Comparing the countmaps on only a subset of qubits")
println("")
display(countmap(samples213))
display(countmap(samples213_full))
println("-----------------------------")
```

    Comparing the full countmaps:
    



    Dict{DitStr{2, 4, BigInt}, Int64} with 4 entries:
      0000 ₍₂₎ => 9276
      1010 ₍₂₎ => 361
      1001 ₍₂₎ => 18
      0011 ₍₂₎ => 345



    Dict{DitStr{2, 4, Int64}, Int64} with 4 entries:
      0000 ₍₂₎ => 9276
      1010 ₍₂₎ => 359
      1001 ₍₂₎ => 9
      0011 ₍₂₎ => 356


    -----------------------------
    
    Comparing the countmaps on only a subset of qubits
    



    Dict{DitStr{2, 3, BigInt}, Int64} with 4 entries:
      000 ₍₂₎ => 9273
      010 ₍₂₎ => 15
      011 ₍₂₎ => 343
      001 ₍₂₎ => 369



    Dict{DitStr{2, 3, Int64}, Int64} with 4 entries:
      000 ₍₂₎ => 9274
      010 ₍₂₎ => 15
      011 ₍₂₎ => 363
      001 ₍₂₎ => 348


    -----------------------------


which looks similar enough to trust that the samples come from the same probability distribution.

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
$$H = J ∑_i^{L-1} X_i X_{i+1} + h ∑_i^L Z_i = U + T.$$
As Ansatz circuits we use the Hamiltonian Variational Ansatz
$$
    U(\vec θ) = ∏_i^p e^{-iθ_{i,U} U} e^{-iθ_{i,T} T} 
$$
with the initial state being the groundstate of the TFIM at $J = 0$, so $|ψ_i⟩ = |0⋯0⟩$


```julia
L = 100 # this is far beyond what is possible with a full wavefunction simulation
J = 1.5
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

ψ_i = FLOYao.zero_state(L);
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

    Iteration 1, energy = -99.94539903107294
    Iteration 2, energy = -100.05435457339475
    Iteration 3, energy = -100.10909645468709
    Iteration 4, energy = -100.17995388555624
    Iteration 5, energy = -100.29171904099248
    Iteration 6, energy = -100.47220887850493
    Iteration 7, energy = -100.76428954038174
    Iteration 8, energy = -101.23517861649259
    Iteration 9, energy = -101.98742997795534
    Iteration 10, energy = -103.1696015011177
    Iteration 11, energy = -104.97820918341799
    Iteration 12, energy = -107.63269450982706
    Iteration 13, energy = -111.29824952351869
    Iteration 14, energy = -115.95044996473993
    Iteration 15, energy = -121.24559138998656
    Iteration 16, energy = -126.53854397808638
    Iteration 17, energy = -131.13147536052566
    Iteration 18, energy = -134.6014550985459
    Iteration 19, energy = -136.93030091440178
    Iteration 20, energy = -138.37167284536923
    Iteration 21, energy = -139.23738100518304
    Iteration 22, energy = -139.7669229470222
    Iteration 23, energy = -140.10252445074727
    Iteration 24, energy = -140.3183531143311
    Iteration 25, energy = -140.45559906861777
    Iteration 26, energy = -140.54217184645134
    Iteration 27, energy = -140.5984013873836
    Iteration 28, energy = -140.6379056962986
    Iteration 29, energy = -140.66876293092122
    Iteration 30, energy = -140.69530934316288
    Iteration 31, energy = -140.71972222959687
    Iteration 32, energy = -140.7430552087384
    Iteration 33, energy = -140.76581206917365
    Iteration 34, energy = -140.78823616608892
    Iteration 35, energy = -140.81044884824954
    Iteration 36, energy = -140.83251393846666
    Iteration 37, energy = -140.8544675259937
    Iteration 38, energy = -140.87633182073378
    Iteration 39, energy = -140.89812172304661
    Iteration 40, energy = -140.91984804138892
    Iteration 41, energy = -140.9415191396324
    Iteration 42, energy = -140.96314182749603
    Iteration 43, energy = -140.9847218716837
    Iteration 44, energy = -141.0062643075656
    Iteration 45, energy = -141.02777364009404
    Iteration 46, energy = -141.04925397961404
    Iteration 47, energy = -141.07070913726895
    Iteration 48, energy = -141.09214269409438
    Iteration 49, energy = -141.11355805226577
    Iteration 50, energy = -141.1349584738442
    Iteration 51, energy = -141.15634711053906
    Iteration 52, energy = -141.17772702688384
    Iteration 53, energy = -141.19910121851333
    Iteration 54, energy = -141.22047262675102
    Iteration 55, energy = -141.2418441503938
    Iteration 56, energy = -141.2632186553526
    Iteration 57, energy = -141.28459898264373
    Iteration 58, energy = -141.30598795510855
    Iteration 59, energy = -141.32738838314748
    Iteration 60, energy = -141.34880306969114
    Iteration 61, energy = -141.3702348145796
    Iteration 62, energy = -141.391686418485
    Iteration 63, energy = -141.41316068647978
    Iteration 64, energy = -141.43466043133293
    Iteration 65, energy = -141.45618847660023
    Iteration 66, energy = -141.4777476595568
    Iteration 67, energy = -141.49934083401257
    Iteration 68, energy = -141.5209708730434
    Iteration 69, energy = -141.54264067166105
    Iteration 70, energy = -141.56435314944287
    Iteration 71, energy = -141.58611125313683
    Iteration 72, energy = -141.60791795925388
    Iteration 73, energy = -141.6297762766571
    Iteration 74, energy = -141.6516892491563
    Iteration 75, energy = -141.67365995811414
    Iteration 76, energy = -141.69569152506713
    Iteration 77, energy = -141.71778711436693
    Iteration 78, energy = -141.73994993584452
    Iteration 79, energy = -141.76218324749823
    Iteration 80, energy = -141.78449035820972
    Iteration 81, energy = -141.80687463048733
    Iteration 82, energy = -141.82933948323847
    Iteration 83, energy = -141.85188839457305
    Iteration 84, energy = -141.87452490463667
    Iteration 85, energy = -141.89725261847698
    Iteration 86, energy = -141.9200752089413
    Iteration 87, energy = -141.94299641960754
    Iteration 88, energy = -141.9660200677507
    Iteration 89, energy = -141.98915004734292
    Iteration 90, energy = -142.0123903320908
    Iteration 91, energy = -142.03574497851153
    Iteration 92, energy = -142.05921812904745
    Iteration 93, energy = -142.0828140152241
    Iteration 94, energy = -142.10653696085234
    Iteration 95, energy = -142.13039138527822
    Iteration 96, energy = -142.15438180668437
    Iteration 97, energy = -142.17851284544648
    Iteration 98, energy = -142.20278922754977
    Iteration 99, energy = -142.22721578807068
    Iteration 100, energy = -142.25179747472944


## Adding support for your own gates

Natively, the only FLO gates that come already shipped with `Yao.jl` are the gates listed [here](#List-of-supported-gates). But there are many more FLO gates, one being for example the `FSWAP` gate which swaps to qubits while making sure to preserve the fermionic commutation relations


```julia
@const_gate FSWAP::ComplexF64 = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 -1]
```

If a gate defines a matrix representation, as we just did for the `FSWAP`gate, `FLOYao` supports them out of the box by manually checking if they are a FLO gate and then computing its matrix representation in the Majorana basis. But this method is fairly slow–though still poly-time and memory–compared to directly implementing `unsafe_apply!(::MajoranaReg, ::YourBlock)` and  `instruct!(::MajoranaReg, ::YourBlock)` and will warn you accordingly


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
    └ @ FLOYao /home/yc20910/PhD/Work/code/FLOYao/src/instruct.jl:56





    MajoranaReg{Float64} with 4 qubits:
    8×8 Matrix{Float64}:
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




    BenchmarkTools.Trial: 6254 samples with 1 evaluation.
     Range (min … max):  713.850 μs …   3.865 ms  ┊ GC (min … max): 0.00% … 74.64%
     Time  (median):     759.798 μs               ┊ GC (median):    0.00%
     Time  (mean ± σ):   794.402 μs ± 234.728 μs  ┊ GC (mean ± σ):  2.27% ±  6.20%
    
       ▃▆▃█▄▄▂                                                       
      ▃███████▇▆▅▄▄▃▃▃▃▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▂▂▂▁▂▁▁▂▂▂▂▂▂▂▂▂▂▂▂ ▃
      714 μs           Histogram: frequency by time         1.24 ms <
    
     Memory estimate: 338.53 KiB, allocs estimate: 495.



To find out what the matrix representation of the `FSWAP` gate in the Majorana basis is, it is easiest to retrace what is happening inside `instruct!(::MajoranaReg, ::AbstractMatrix, locs)`
You can use


```julia
@which instruct!(mreg, mat(FSWAP), (1,2))
```




instruct!(reg::<b>MajoranaReg</b>, gate::<b>AbstractMatrix</b>, locs) in FLOYao at <a href="https://github.com/PhaseCraft/FLOYao.jl/tree/4d66a3d6b5629ae5eb391e1a1b585c7abd1a1d53//src/instruct.jl#L49" target="_blank">/home/yc20910/PhD/Work/code/FLOYao/src/instruct.jl:49</a>



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




    BenchmarkTools.Trial: 10000 samples with 585 evaluations.
     Range (min … max):  188.162 ns …   4.742 μs  ┊ GC (min … max): 0.00% … 92.38%
     Time  (median):     213.903 ns               ┊ GC (median):    0.00%
     Time  (mean ± σ):   249.819 ns ± 296.990 ns  ┊ GC (mean ± σ):  9.48% ±  7.55%
    
         ▅▄█▆▅▆▆▅▄▄▄▃▂▁▁                               ▁▁▁▁         ▂
      ▃▅█████████████████▇▆▅▆▄▄▄▅▃▄▅▅▅▃▄▅▅▅▆▆▆▆▆▇▇▇▇████████████▇▇▆ █
      188 ns        Histogram: log(frequency) by time        388 ns <
    
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
```math
\begin{aligned}
    ⟨ψ|UHU^†|ψ⟩ &= \frac{i}{4} \tilde H^{mn} ⟨Ω|γ_{m} γ_{n}|Ω⟩ \\
                &= - \frac{1}{2} ∑_{i} \tilde H^{2i-1,2i} \\
                &= - \frac{1}{2} ∑_{i>k} R^{2i-1}_{m} R^{2i}_{n} H^{mn} \\
\end{aligned}.
```
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

