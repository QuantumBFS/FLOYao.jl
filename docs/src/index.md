# FLOYao


```@docs
FLOYao.FLOYao
```

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

and define a (here similarly arbitrary) observable

```jldoctest quickstart; output=false
hamiltonian = (kron(nq, 1=>Z) + 2kron(nq, 1=>X, 2=>Z, 3=>Z, 4=>X) + 3.5put(nq, 2=>Z)
               + 0.5kron(nq, 1=>Z, 2=>Z) - kron(nq, 2 => X, 4 => Y))

# output
nqubits: 4
+
├─ +
│  ├─ kron
│  │  └─ 1=>Z
│  ├─ [scale: 2] kron
│  │     ├─ 1=>X
│  │     ├─ 2=>Z
│  │     ├─ 3=>Z
│  │     └─ 4=>X
│  ├─ [scale: 3.5] put on (2)
│  │     └─ Z
│  └─ [scale: 0.5] kron
│        ├─ 1=>Z
│        └─ 2=>Z
└─ [-] kron
      ├─ 2=>X
      └─ 4=>Y
```

and finally create a register in the computational zero state via

```jldoctest quickstart
reg = FLOYao.zero_state(nq)

# output
MajoranaReg{Float64} with 4 qubits:
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
 -1.0  -0.0                 -0.0                 …  -0.0                 -0.0
 -0.0  -0.9238795325112867   0.3826834323650898     -0.0                 -0.0
  0.0   0.3826834323650898   0.9238795325112867      0.0                  0.0
  0.0   0.0                  0.0                     0.3826834323650898   0.0
  0.0   0.0                  0.0                     0.0                  0.0
  0.0   0.0                  0.0                 …   0.0                  0.0
  0.0   0.0                  0.0                     0.9238795325112867   0.0
  0.0   0.0                  0.0                     0.0                  1.0
```

and the same goes for expectation values of observables

```jldoctest quickstart
expect(hamiltonian, reg => circuit)

# output
5.04883281702193
```

or even gradients of these expectation values with respect to the circuit parameters

```jldoctest quickstart
state_grad, params_grad = expect'(hamiltonian, reg => circuit)

# output
MajoranaReg{Float64}(4) => [0.0, -0.9130135182550002, 0.20220773043883133, -5.551115123125783e-17]
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


