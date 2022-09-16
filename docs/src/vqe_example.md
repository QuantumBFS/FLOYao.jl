# Example: VQE for the transverse field Ising model

One prime example of a free fermion Hamiltonian and FLO circuit is 
the variational quantum eigensolver with the hamiltonian variational ansatz 
applied  the transverse field Ising model. It is a good example demonstrating
the capabilities of `FLOYao`.

The Hamiltonian is given as 
```math
    H = J ∑_i^{L-1} X_i X_{i+1} + h ∑_i^L Z_i = U + T.
```
and as Ansatz circuits we use the Hamiltonian Variational Ansatz
```math
    U(\vec θ) = ∏_i^p e^{-iθ_{i,U} U} e^{-iθ_{i,T} T} 
```
with the initial state being the groundstate of the TFIM at $J = 0$, so $|ψ_i⟩ = |0 ⋯ 0⟩$.

First, we define the Hamiltonian

```jldoctest tfimvqe
using FLOYao, Yao
L = 100 # this is far beyond what is possible with a full wavefunction simulation
J = 1.5 
h = -1.
p = 10  # number of VQE layers
U = map(1:L-1) do i
    J * kron(L, i => X, i+1 => X)
end |> sum

T = map(1:L) do i
    h * kron(L, i => Z)
end |> sum

hamiltonian = T + U
# not really needed, but here to circumvent some doctest  restrictions
summary(hamiltonian)

# output
"Add{2}"
```


and the ansatz circuit

```jldoctest tfimvqe; output=false
circuit = chain(L)
for _ in 1:p
    for i in 1:L-1
        push!(circuit, rot(kron(L, i => X, i+1 => X), 0.))
    end
    for i in 1:L
        push!(circuit, put(L, i => Rz(0.)))
    end
end

# output
```

as well as the initial state

```jldoctest tfimvqe
reg = FLOYao.zero_state(L)
typeof(reg)

# output
MajoranaReg{Float64}
```

Now that we defined the hamiltonian, the ansatz circuit and the initial state
we can perform simple gradient descent on the energy expectation value to find
an approximation to the groundstate of $H$:


```jldoctest tfimvqe
iterations = 100
gamma = 2e-2

# set the initial parameters
nparams = nparameters(circuit)
dispatch!(circuit, ones(nparams) ./ 100) # fix initial parameters for reproducibility

for i in 1:iterations
    _, grad = expect'(hamiltonian, reg => circuit)
    dispatch!(-, circuit, gamma * grad)
    println("Iteration $i, energy = $(round(expect(hamiltonian, reg => circuit), digits=4))")
end

# output
Iteration 1, energy = -99.7577
Iteration 2, energy = -100.1858
Iteration 3, energy = -100.3843
Iteration 4, energy = -100.6364
Iteration 5, energy = -101.0342
Iteration 6, energy = -101.6735
Iteration 7, energy = -102.6934
Iteration 8, energy = -104.2918
Iteration 9, energy = -106.7225
Iteration 10, energy = -110.2443
Iteration 11, energy = -114.983
Iteration 12, energy = -120.7112
Iteration 13, energy = -126.7071
Iteration 14, energy = -131.975
Iteration 15, energy = -135.8116
Iteration 16, energy = -138.1614
Iteration 17, energy = -139.4181
Iteration 18, energy = -140.0345
Iteration 19, energy = -140.3264
Iteration 20, energy = -140.4677
Iteration 21, energy = -140.5422
Iteration 22, energy = -140.5877
Iteration 23, energy = -140.6206
Iteration 24, energy = -140.648
Iteration 25, energy = -140.6729
Iteration 26, energy = -140.6966
Iteration 27, energy = -140.7196
Iteration 28, energy = -140.7423
Iteration 29, energy = -140.7647
Iteration 30, energy = -140.7869
Iteration 31, energy = -140.809
Iteration 32, energy = -140.831
Iteration 33, energy = -140.8529
Iteration 34, energy = -140.8747
Iteration 35, energy = -140.8965
Iteration 36, energy = -140.9182
Iteration 37, energy = -140.9398
Iteration 38, energy = -140.9614
Iteration 39, energy = -140.983
Iteration 40, energy = -141.0045
Iteration 41, energy = -141.026
Iteration 42, energy = -141.0475
Iteration 43, energy = -141.0689
Iteration 44, energy = -141.0904
Iteration 45, energy = -141.1118
Iteration 46, energy = -141.1332
Iteration 47, energy = -141.1545
Iteration 48, energy = -141.1759
Iteration 49, energy = -141.1973
Iteration 50, energy = -141.2187
Iteration 51, energy = -141.24
Iteration 52, energy = -141.2614
Iteration 53, energy = -141.2828
Iteration 54, energy = -141.3042
Iteration 55, energy = -141.3256
Iteration 56, energy = -141.347
Iteration 57, energy = -141.3685
Iteration 58, energy = -141.3899
Iteration 59, energy = -141.4114
Iteration 60, energy = -141.4329
Iteration 61, energy = -141.4544
Iteration 62, energy = -141.476
Iteration 63, energy = -141.4976
Iteration 64, energy = -141.5193
Iteration 65, energy = -141.5409
Iteration 66, energy = -141.5627
Iteration 67, energy = -141.5844
Iteration 68, energy = -141.6062
Iteration 69, energy = -141.6281
Iteration 70, energy = -141.65
Iteration 71, energy = -141.672
Iteration 72, energy = -141.694
Iteration 73, energy = -141.7161
Iteration 74, energy = -141.7383
Iteration 75, energy = -141.7605
Iteration 76, energy = -141.7829
Iteration 77, energy = -141.8052
Iteration 78, energy = -141.8277
Iteration 79, energy = -141.8503
Iteration 80, energy = -141.8729
Iteration 81, energy = -141.8956
Iteration 82, energy = -141.9184
Iteration 83, energy = -141.9414
Iteration 84, energy = -141.9644
Iteration 85, energy = -141.9875
Iteration 86, energy = -142.0107
Iteration 87, energy = -142.0341
Iteration 88, energy = -142.0576
Iteration 89, energy = -142.0812
Iteration 90, energy = -142.1049
Iteration 91, energy = -142.1287
Iteration 92, energy = -142.1527
Iteration 93, energy = -142.1768
Iteration 94, energy = -142.2011
Iteration 95, energy = -142.2255
Iteration 96, energy = -142.25
Iteration 97, energy = -142.2748
Iteration 98, energy = -142.2997
Iteration 99, energy = -142.3247
Iteration 100, energy = -142.3499
```

Hopefully, this is a good enough approximation to the groundstate. We can now 
use this to sample from the state we found:

```jldoctest tfimvqe
using Random
samples = measure(reg |> circuit, nshots=10, rng=MersenneTwister(42))

# output

10-element Vector{DitStr{2, 100, BigInt}}:
 0000000000000000000000000000000000000000000000000011000000000000000000000101000000000000111100110000 ₍₂₎
 0000000000000000000000000000000000000000000000000000000000000000000000000011011001100000000000000000 ₍₂₎
 0000000000000000000000000000000000000000000000011000000000001100000000110000000000001100000110000000 ₍₂₎
 0000000000000000000000000000000000000011000000011000000000000110000000000011000000000000000011000011 ₍₂₎
 0000000000000000000000000000000000000000000000000000000000000000000001100000000000011111100000000000 ₍₂₎
 0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000011000011 ₍₂₎
 0000000000000000000000000000000000000000000000001010000000101000000000000000000000011000000000000000 ₍₂₎
 0000000000000000000000000000000000000000000000000000000001100000000000000000000000000001100000000000 ₍₂₎
 0000000000000000000000000000000000000000000000000000110000000010100000000001100000000000000011110000 ₍₂₎
 0000000000000000000000000000000000000110000000000010001000000000110000000000001100000011011000000000 ₍₂₎
```
