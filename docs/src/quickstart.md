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
    println("Iteration $i, energy = $(expect(hamiltonian, reg => circuit))")
end

# output

Iteration 1, energy = -99.75768874435411
Iteration 2, energy = -100.18576067696775
Iteration 3, energy = -100.38432173792287
Iteration 4, energy = -100.63643044122702
Iteration 5, energy = -101.03423092444348
Iteration 6, energy = -101.67354500209733
Iteration 7, energy = -102.69338647473786
Iteration 8, energy = -104.29178300816409
Iteration 9, energy = -106.72246768352733
Iteration 10, energy = -110.24427073980664
Iteration 11, energy = -114.9829856652035
Iteration 12, energy = -120.7112264427376
Iteration 13, energy = -126.70705586127221
Iteration 14, energy = -131.97496218581824
Iteration 15, energy = -135.811588520062
Iteration 16, energy = -138.16137020312354
Iteration 17, energy = -139.418101921114
Iteration 18, energy = -140.0344906389265
Iteration 19, energy = -140.32644601982946
Iteration 20, energy = -140.46768708932973
Iteration 21, energy = -140.54215336725701
Iteration 22, energy = -140.58765045468903
Iteration 23, energy = -140.62059971371164
Iteration 24, energy = -140.64802527737191
Iteration 25, energy = -140.67292942947097
Iteration 26, energy = -140.6966109232686
Iteration 27, energy = -140.71964633900248
Iteration 28, energy = -140.74230283171434
Iteration 29, energy = -140.76471213063707
Iteration 30, energy = -140.78694464194248
Iteration 31, energy = -140.80904160883958
Iteration 32, energy = -140.83102948670762
Iteration 33, energy = -140.85292664660233
Iteration 34, energy = -140.87474668214685
Iteration 35, energy = -140.8965001557173
Iteration 36, energy = -140.91819559107984
Iteration 37, energy = -140.93984007923905
Iteration 38, energy = -140.96143967164957
Iteration 39, energy = -140.98299964815365
Iteration 40, energy = -141.00452470633957
Iteration 41, energy = -141.02601909898442
Iteration 42, energy = -141.04748673578447
Iteration 43, energy = -141.06893125978692
Iteration 44, energy = -141.09035610549105
Iteration 45, energy = -141.11176454346037
Iteration 46, energy = -141.13315971488998
Iteration 47, energy = -141.1545446586167
Iteration 48, energy = -141.17592233241598
Iteration 49, energy = -141.1972956299425
Iteration 50, energy = -141.21866739434242
Iteration 51, energy = -141.2400404293133
Iteration 52, energy = -141.26141750819335
Iteration 53, energy = -141.28280138153886
Iteration 54, energy = -141.30419478353016
Iteration 55, energy = -141.3256004374722
Iteration 56, energy = -141.34702106059856
Iteration 57, energy = -141.36845936833754
Iteration 58, energy = -141.38991807816123
Iteration 59, energy = -141.41139991312326
Iteration 60, energy = -141.4329076051475
Iteration 61, energy = -141.4544438981408
Iteration 62, energy = -141.4760115509648
Iteration 63, energy = -141.4976133403108
Iteration 64, energy = -141.51925206349878
Iteration 65, energy = -141.54093054123172
Iteration 66, energy = -141.56265162031173
Iteration 67, energy = -141.58441817634179
Iteration 68, energy = -141.60623311641493
Iteration 69, energy = -141.62809938180553
Iteration 70, energy = -141.65001995066322
Iteration 71, energy = -141.67199784071812
Iteration 72, energy = -141.69403611199868
Iteration 73, energy = -141.71613786956488
Iteration 74, energy = -141.73830626625946
Iteration 75, energy = -141.76054450547576
Iteration 76, energy = -141.78285584394717
Iteration 77, energy = -141.80524359455376
Iteration 78, energy = -141.82771112915003
Iteration 79, energy = -141.8502618814116
Iteration 80, energy = -141.87289934970264
Iteration 81, energy = -141.8956270999624
Iteration 82, energy = -141.91844876861146
Iteration 83, energy = -141.94136806548082
Iteration 84, energy = -141.96438877675814
Iteration 85, energy = -141.98751476795846
Iteration 86, energy = -142.0107499869182
Iteration 87, energy = -142.03409846681043
Iteration 88, energy = -142.05756432919154
Iteration 89, energy = -142.08115178707055
Iteration 90, energy = -142.10486514801698
Iteration 91, energy = -142.12870881729827
Iteration 92, energy = -142.1526873010611
Iteration 93, energy = -142.17680520955616
Iteration 94, energy = -142.20106726040987
Iteration 95, energy = -142.22547828195925
Iteration 96, energy = -142.25004321664372
Iteration 97, energy = -142.2747671244732
Iteration 98, energy = -142.299655186575
Iteration 99, energy = -142.32471270883173
Iteration 100, energy = -142.3499451256187
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
