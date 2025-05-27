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

# ------------------------------
# Calculating expectation values
# ------------------------------
function Yao.expect(ham::MajoranaSum{Complex{HT}}, reg::MajoranaReg{RT}) where {HT, RT}
    C = covariance_matrix(reg)
    even_terms = Iterators.filter(iseven ∘ length, ham)
    sum(even_terms, init=zero(promote_type(HT, RT))) do term
        l = length(term) ÷ 2
        coeff = real(1im^l * term.coeff)
        l == 0 ? coeff : coeff * pfaffian!(C[term.indices, term.indices])
    end
end

"""
    expect(op::Union{AbstractBlock,MajoranaSum}, reg::MajoranaReg)
    expect(op::Union{AbstractBlock,MajoranaSum}, (reg => circuit)::Pair{<:MajoranaReg,<:AbstractBlock})

Get the expectation value of an operator, the second parameter can be a
register `reg` or a pair of input register and circuit `reg => circuit`.

    expect'(op::Union{AbstractBlock,MajoranaSum}, reg => circuit::) -> Pair{<:MajoranaReg,<:AbstractVector}
    expect'(op::Union{AbstractBlock,MajoranaSum}, reg::MajoranaReg) -> MajoranaReg

Obtain the gradient with respect to registers and circuit parameters. For pair
input, the second return value is a pair of `gψ => gparams`, with `gψ` the gradient
of input state and `gparams` the gradients of circuit parameters. For register
input, the return value is a register.
"""
function Yao.expect(op::AbstractBlock, reg::MajoranaReg)
    YaoBlocks._check_size(reg, op)
    majoranasum = yaoblock2majoranasum(eltype(reg), op)
    return expect(majoranasum, reg)
end

# ----------------------
# Calculating fidelities
# ----------------------
"""
    bitstring_probability(reg::MajoranaReg, bit_string::BitStr)

The probability to measure a given `bit_string` when measuring `reg`

# Todo
 - Rewrite this using pfaffians
"""
function bitstring_probability(reg::MajoranaReg{T}, bit_string::DitStr{2,N,ST}) where {T,N,ST}
    @assert nqubits(reg) == N
    p = one(T)
    M = covariance_matrix(reg)
    for i in 1:N
        ni = bit_string[i]
        pi = (1 + (-1)^ni * M[2i-1, 2i]) / 2
        pi ≈ 0 && return zero(T)
        p *= pi
        update_covariance_matrix!(M, i, pi, ni)
    end
    # floating point issues can cause very small probabilities to get negative.
    return p > zero(T) ? p : zero(T)
end

"""
    fidelity(reg1::MajoranaReg, reg2::MajoranaReg)
    fidelity'(pair_or_majoranareg, pair_or_majoranareg) -> (g1, g2)

The fidelity ``|⟨ψ|φ⟩|`` between the FLO states ``|ψ⟩`` and ``|φ⟩`` defined
by `reg1` and `reg2`.

The gradient with respect to registers and circuit parameters. For a pair input
`reg => circuit` the returned gradient is a pair of `greg => gparams` with
`greg` the gradient with respect to the input state and `gparams` the gradient
with respect to circuit parameters. For a `reg` input the return value is the
gradient with respect to the register.

!!! note
    This definition differs from the standard definition ``|⟨φ|ψ⟩|²`` by a square 
    root, but is consistent with the one used in `Yao.jl`
"""
function Yao.fidelity(reg1::MajoranaReg{T1}, reg2::MajoranaReg{T2}) where {T1,T2}
    T = promote_type(T1, T2)
    nq = nqubits(reg1)
    reg = MajoranaReg(reg2.state' * reg1.state)
    return √bitstring_probability(reg, BitStr{nq}(0))
end
