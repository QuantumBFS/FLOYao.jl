# [Mathematical Background](@id background)

This section is here, more to fix the convention of
[Jordan-Wigner transform](https://en.wikipedia.org/wiki/Jordan%E2%80%93Wigner_transformation)
and [Majorana operators](https://en.wikipedia.org/wiki/Majorana_fermion) that we use here,
and less to explain the full theory behind those. For the latter, we, once again, recommend [Classical simulation of noninteracting-fermion quantum circuits](https://arxiv.org/abs/quant-ph/0108010).

We define the Majorana operators $γ_i$ via 
```math
    γ_{2i-1} = ∏_{j=1}^{i-1} (-Z_j) X_i
    \qquad \textrm{and} \qquad
    γ_{2i} = -∏_{j=1}^{i-1} (-Z_j) Y_i.
```
This implies the normal fermionic creation and annihilation operators are given by
```math
    c_j = \frac{1}{2} (γ_{2j-1} + iγ_{2j})
    \quad \textrm{and} \quad
    c_j^† = \frac{1}{2} (γ_{2j-1} - iγ_{2j})
```
and products of two Majorana operators are of the form
```math
    σ_i \left(∏_{i<j<k} -Z_j \right) σ_k
    \quad \textrm{or} \quad
    Z_i
```
with $σ_i, σ_k ∈ \{X, Y\}$.

Any unitary that takes all Majorana operators to a linear combination of
Majorana operators under conjugation, i.e. that satisfies
```math
    U γ_i U^† = R_i^j γ_j
```
with some $R ∈ O(2n)$ is a FLO unitary. In particular, if a unitary is of the form 
```math
    U = e^{-iθH}
```
with 
```math
    H = \frac{i}{4} \sum_{i,j} H^{ij} γ_i γ_j
```
it is a FLO unitary with $R ∈ SO(2n)$.

But note, that not all FLO unitaries are of that form. For example, $X_i$ is
also a FLO gate since it either commutes or anti-commutes with all Majorana
operators, but the associated matrix $R$ always has determinant $-1$.

Calculating the expectation values of hamiltonians like the one above when
evolving the vacuum state with FLO circuits is efficiently possible. First
evolve the Hamiltonian in the Heisenber picture to
```math
    UHU^† = \frac{i}{4} R^{m}_{i} R^{n}_{j} H^{ij} γ_{m} γ_{n} 
           =: \frac{i}{4} \tilde H^{mn} γ_{m} γ_{n}.
```
and then compute the expectation value
```math
\begin{aligned}
    ⟨ψ|UHU^†|ψ⟩ &= \frac{i}{4} \tilde H^{mn} ⟨Ω|γ_{m} γ_{n}|Ω⟩ \\
                &= - \frac{1}{2} ∑_{i} \tilde H^{2i-1,2i} \\
                &= - \frac{1}{2} ∑_{i} R^{2i-1}_{m} R^{2i}_{n} H^{mn} \\
\end{aligned}.
```
From the first to second line one needs to carefully think which of the 
$⟨Ω|γ_{m} γ_{n}|Ω⟩$ are zero and which cancel each other out due to the anti-symmetry of $H^{mn}$.

