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
    M = zeros(eltype(reg), 2nq, 2nq)
    G = I(nq) ⊗ [0 1; -1 0]
    return reg.state * G * reg.state'
end

"""
Updates the covariance matrix `M`, given that the `i`th qubit was measured
in the state `ni` with probability `pi`.

# Note: After this procedure only the entries M[p,q] with 2i-1 < p < q
will be correct. This is sufficient to calculate measurement probabilities
for j > i.
"""
function update_covariance_matrix!(M, i, pi, ni)
    n = size(M,1)
    for p in 2i+1:n 
        for q in p+1:n
            M[p,q] += (-1)^ni * M[2i-1,q] * M[2i,p] / (2pi)
            M[p,q] -= (-1)^ni * M[2i-1,p] * M[2i,q] / (2pi)
        end
    end
    return M
end

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
@inline function sample(covmat, locs, ids,
                        rng=Random.GLOBAL_RNG)
    out = BigInt(0)
    nq = length(locs)
    for i in 1:nq
        pi = (1 + covmat[2locs[ids[i]]-1,2locs[ids[i]]]) / 2
        ni = rand(rng) > pi
        out += ni * 2^(ids[i]-1)
        update_covariance_matrix!(covmat, i, ni ? 1-pi : pi, ni)
    end
    return out
end

# a slightly faster version, when all qubits get sampled in their normal
# order
@inline function sample(covmat, rng=Random.GLOBAL_RNG)
    out = BigInt(0)
    nq = size(covmat,1) ÷ 2
    for i in 1:nq
        pi = (1 + covmat[2i-1,2i]) / 2
        ni = rand(rng) > pi
        out += ni * 2^(i-1)
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
function Yao.measure(reg::MajoranaReg, locs;
                     nshots::Int=1, rng=Random.GLOBAL_RNG)
    covmat = covariance_matrix(reg)
    ids = sortperm(locs)
    samples = Vector{BitStr{length(locs),BigInt}}(undef, nshots)
    M = similar(covmat)
    for i in 1:nshots
        M .= covmat
        samples[i] = sample(M, locs, ids, rng)
    end
    return samples
end

function Yao.measure(reg::MajoranaReg;
                     nshots::Int=1, rng=Random.GLOBAL_RNG)
    covmat = covariance_matrix(reg)
    samples = Vector{BitStr{size(covmat,1)÷2,BigInt}}(undef, nshots)
    M = similar(covmat)
    for i in 1:nshots
        M .= covmat
        samples[i] = sample(M, rng)
    end
    return samples
end

"""
    measure([postprocess::Union{NoPostProcess, ResetTo},] reg::MajoranaReg; rng=Random.GLOBAL_RNG)

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
function Yao.measure!(reg::MajoranaReg; rng=Random.GLOBAL_RNG)
    nq = nqubits(reg)
    covmat = covariance_matrix(reg)
    res = BitStr{nq}(sample(covmat, rng))
    FLOYao.product_state!(reg, res)
    return res
end

function Yao.measure!(postprocess::NoPostProcess, reg::MajoranaReg; rng=Random.GLOBAL_RNG)
    res = measure!(reg, rng=rng)
    return res
end

function Yao.measure!(postprocess::ResetTo, reg::MajoranaReg; rng=Random.GLOBAL_RNG)
    nq = nqubits(reg)
    covmat = covariance_matrix(reg)
    res = BitStr{nq}(sample(covmat, rng))
    FLOYao.product_state!(reg, postprocess.x)
    return res
end
