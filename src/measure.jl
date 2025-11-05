#=
#  Authors:   Jan Lukas Bosse
#  Copyright: 2022 Phasecraft Ltd.
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
=#

# ------------------------------------------------
# Measuring registers and sampling from FLO states
# ------------------------------------------------
"""
    covariance_matrix(reg::MajoranaReg)

The covariance matrix
```math
    M_{pq} = \\frac{i}{2} ⟨Ω|U^†[γ_p,γ_q]U|Ω⟩
```
of a FLO state ``U|Ω⟩``.
"""
function covariance_matrix(reg::MajoranaReg)
    nq = nqubits(reg)
    # TODO: Don't instantiate G, but implement R * G by row (col?) swapping
    G = kron(I(nq), [1im 1; -1 1im])
    M = reg.state * G * reg.state'
    M = (M - transpose(M)) / 2
    M[diagind(M)] .= 1im
    return M
end

"""
Updates the covariance matrix `M`, given that the `i`th qubit was measured
in the state `ni` with probability `pi`.
"""
function update_covariance_matrix!(M::AbstractMatrix{T}, i, pi, ni) where {T}
    n = size(M, 1)
    m = copy(M)
    for p in 1:n, q in 1:n
        if ((p + 1) ÷ 2 == i) ⊻ ((q + 1) ÷ 2 == i)
            M[p,q] = zero(T)
        else
            M[p,q] -= (-1)^ni * m[2i-1,p] * m[2i,q] / (2pi)
            M[p,q] += (-1)^ni * m[2i-1,q] * m[2i,p] / (2pi)
        end
    end

    # manually fix the 2i-2, 2i entries. No idea why this is neccessary...
    # It shouldn't be?
    # M[2i-1, 2i] = (-1)^ni
    # M[2i, 2i-1] = -(-1)^ni

    return M
end

"""
    state_matrix(cov_mat::AbstractMatrix) -> AbstractMatrix

Construct a state matrix ∈ O(n) from a covariance_matrix.
This is inverse to [`covariance_matrix`](@ref).
"""
function state_matrix(cov_mat::AbstractMatrix{T}) where {T}
    nq = size(cov_mat, 1) ÷ 2
    cov_mat[diagind(cov_mat)] .= zero(T)
    q, h = hessenberg(cov_mat)
    state = real.(Matrix(q))
    for i in 1:nq
        if real(h[2i, 2i-1]) > zero(real(T))
            state[:, 2i] .*= -1
        end
    end
    return state
end


# TODO: Is `ids` actually needed with the new `update_covariance_matrix`?
"""
    sample!(covmat, locs=1:size(covmat,1)÷2, ids=sortperm(locs), rng=GLOBAL_RNG)

Take a computational basis state sample from a FLO state with given `cov`arianve` mat`rix

# Arguments
 - `covmat`: The 2n×2n covariance matrix ``M_{ij} = i ⟨ψ|γ_i γ_j|ψ⟩ - i δ_{ij}``.
             This matrix will be overwritten by the function!
 - `locs=1:size(covmat,1)÷2`: The list of qubits to measure
 - `ids=sortperm(locs)`: The sorting permutation for `locs`. You don't need 
                         to pass this, but precomputing it can make this faster
 - `rng=Random.GLOBAL_RNG`: The random number generator to use for sampling.

# Warning
This function works inplace on the covmat argument to avoid allocations.
If you need multiple samples from the same covmat, make sure to pass a copy in.
"""
@inline function sample!(covmat, locs, rng=Random.GLOBAL_RNG)
    out = BigInt(0)
    nq = length(locs)
    display(round.(covmat, digits=2))
    for i in locs
        pi = (1 + covmat[2i-1,2i]) / 2
        @show pi
        pi = real(pi)
        ni = rand(rng) > pi
        out += ni * 2^(i-1)
        update_covariance_matrix!(covmat, i, ni ? 1-pi : pi, ni)
        display(round.(covmat, digits=2))
    end
    return out
end

# a slightly faster version, when all qubits get sampled in their normal
# order
# TODO: Delete this when I don't need it anymore
@inline function sample!(covmat, rng=Random.GLOBAL_RNG)
    out = BigInt(0)
    nq = size(covmat,1) ÷ 2
    for i in 1:nq
        pi = (1 + covmat[2i-1,2i]) / 2
        ni = rand(rng) > pi
        out += ni * BigInt(2)^(i-1)
        update_covariance_matrix!(covmat, i, ni ? 1-pi : pi, ni)
    end
    return out
end

"""
    measure(reg::MajoranaReg[, locs]; nshots=1, rng=Random.GLOBAL_RNG) -> Vector{BitStr}

Measure a MajoranaReg and return the results as a vector of `BitStr`ings. This
is a cheating version of `measure!` that does not need to recompute the `reg` 
for each new sample and also does not alter it.

# Arguments
 - `reg::MajoranaReg`: The state register
 - `locs`: The list of qubits to measure. Defaults to all qubits.
 - `nshots=1`: The number of samples to take
 - `rng=Random.GLOBAL_RNG`: The random number generator to use
"""
function Yao.measure(reg::MajoranaReg;
                     nshots::Int=1, rng=Random.GLOBAL_RNG)
    covmat = covariance_matrix(reg)
    samples = Vector{BitStr{size(covmat,1)÷2,BigInt}}(undef, nshots)
    M = similar(covmat)
    for i in 1:nshots
        M .= covmat
        samples[i] = sample!(M, rng)
    end
    return samples
end

function Yao.measure(reg::MajoranaReg, locs;
                     nshots::Int=1, rng=Random.GLOBAL_RNG)
    covmat = covariance_matrix(reg)
    ids = sortperm(locs)
    samples = Vector{BitStr{length(locs),BigInt}}(undef, nshots)
    M = similar(covmat)
    for i in 1:nshots
        M .= covmat
        samples[i] = sample!(M, locs, ids, rng)
    end
    return samples
end

"""
    measure!([postprocess,] reg::MajoranaReg[, locs]; rng=Random.GLOBAL_RNG)

Measure a Majorana `reg`ister in the computational basis and return the 
resulting `BitStr`.

# Arguments
 - `postprocess`: Is the post processing method
    - `NoPostProcess()` (the default) will collapse the state to the measurement
                        outcome state.
    - `ResetTo(bit_str)` will reset the state to the product state specified by 
                         `bit_str`
 - `reg::MajoranaReg`: The quantum register
 - `rng::Random.GLOBAL_RNG`: The RNG to use
"""
function Yao.measure!(reg::MajoranaReg, locs=1:nqubits(reg); rng=Random.GLOBAL_RNG)
    nq = length(locs)
    @show det(reg.state)
    covmat = covariance_matrix(reg)
    res = BitStr{nq}(sample!(covmat, locs, rng))
    display(round.(covmat, digits=2))
    reg.state .= state_matrix(covmat)
    return res
end

function Yao.measure!(postprocess::NoPostProcess, reg::MajoranaReg, locs=1:nqubits(reg);
                      rng=Random.GLOBAL_RNG)
    return measure!(reg, locs, rng=rng)
end

function Yao.measure!(postprocess::ResetTo, reg::MajoranaReg, locs=1:nqubits(reg);
                      rng=Random.GLOBAL_RNG)
    res = measure!(reg, locs, rng=rng)
    FLOYao.product_state!(reg, postprocess.x)
    return res
end
