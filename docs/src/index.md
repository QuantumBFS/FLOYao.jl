# FLOYao

A [Yao.jl](https://github.com/QuantumBFS/Yao.jl) backend to efficiently simulated
fermionic linear optics (FLO) circuits in  based on
[Classical simulation of noninteracting-fermion quantum circuits](https://arxiv.org/abs/quant-ph/0108010)
and [Disorder-assisted error correction in Majorana chains](https://arxiv.org/abs/1108.3845).
FLO circuits are a class of quantum circuits that are closely related to
non-interacting fermions and can be efficiently simulated on classical
computers, similar to the way Clifford circuits can be efficiently classically
simulated, as is done in [YaoClifford.jl](https://github.com/QuantumBFS/YaoClifford.jl).

A quick introduction to fermionic linear optics circuits is found in
[Background](@ref background) section and a more in-depth introduction in e.g. the two
papers linked above.

## Installation
`FLOYao` can be simply installed from the REPL via

```jl-repl
pkg> add FLOYao
```

## Quickstart
First import `FLOYao` and `Yao`
```jldoctest quickstart; output=false
using FLOYao, Yao

# output
```

then build a (here somewhat arbitrary) circuit consisting only of [Supported gates](@ref)

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

circuit

# output
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
```

and define an observable that is a sum of squares of [Majorana operators](@ref background)

```jldoctest quickstart; output=false
hamiltonian = xxg1 + xxg2 + xxg3 + kron(nq, 2=>Z) + kron(nq, 3=>Z)

# output
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
```

and finally create a register in the computational zero state via

```jldoctest quickstart
reg = FLOYao.zero_state(nq)

# output
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
```

Applying the circuit to the register works then exactly the same way as for a normal `ArrayReg` register:

```jldoctest quickstart
apply(reg, circuit)

# output
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
```

and the same goes for expectation values of observables

```jldoctest quickstart
expect(hamiltonian, reg => circuit)

# output
1.8535533905932737
```

or even gradients of these expectation values with respect to the circuit parameters

```jldoctest quickstart
state_grad, params_grad = expect'(hamiltonian, reg => circuit)

# output
MajoranaReg{Float64}(4) => [0.0, -0.3535533905932738, -0.3535533905932738, 0.0]
```

## Contents

```@contents
Pages = [
    "vqe_example.md",
    "features/features.md",
    "features/supported_gates.md",
    "background.md",
    "adding_gates.md",
    "known_restrictions.md"
]
```


