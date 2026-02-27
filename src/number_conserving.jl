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
function check_definite_occupation(J::AbstractMatrix{T}, atol=1e-10) where T
    n = size(J, 1) ÷ 2
    return all(Iterators.product(1:n, 1:n)) do (i, j)
        return if i != j
            (isapproxzero(J[2i-1, 2j-1], atol=atol) && isapproxzero(J[2i-1, 2j], atol=atol)
             && isapproxzero(J[2i, 2j-1], atol=atol) && isapproxzero(J[2i, 2j], atol=atol))
        else
            (isapproxzero(J[2i-1, 2j-1], atol=atol) && isapprox(abs(J[2i-1, 2j]), one(T))
             && isapprox(abs(J[2i, 2j-1]), one(T)) && isapproxzero(J[2i, 2j], atol=atol))
        end
    end
end

"""
    initially_occupied_sites(reg::MajoranaReg)

Find the qubits on which an X gate was applied to in the FLO state `reg` with
definite particle number.
"""
function initially_occupied_sites(reg::MajoranaReg)
    nq = nqudits(reg)
    J = FLOYao.covariance_matrix(FLOYao.zero_state(nq))
    R = state(reg)
    J_prime = R' * J * R
    check_definite_occupation(J_prime) || throw(IndefiniteOccupationException("reg does not have a definite particle number"))
    return filter(1:nq) do i
        J_prime[2i-1, 2i] < 0
    end
end

"""
    xlocs_and_number_conserving_unitary(reg::MajoranaReg)

Decompose the unitary that prepared `reg` from the zero state
into the `xlocs` and the number conserving unitary.
"""
function productstate_and_circuit_decomposition(reg::MajoranaReg)
    nq = nqubits(reg)
    locs = initially_occupied_sites(reg)
    x_state = FLOYao.zero_state_like(reg)
    x_locs = initially_occupied_sites(reg) 
    foreach(x_locs) do qb
        apply!(x_state, put(nq, qb => X))
    end
    return x_locs, state(reg) * state(x_state)'
end

"""
    number_conserving_unitary(O::AbstractMatrix)

Convert the orthogonal `2n×2n` matrix `O` into the equivalent `n×n` unitary `U`.
Note that this only works of `O` corresponds to a number conserving FLO circuit.
"""
function majorana2unitary(O::AbstractMatrix{T}) where {T<:Real}
    n = size(O, 1) ÷ 2
    U = zeros(complex(T), n, n)
    for (i, j) in Iterators.product(1:n, 1:n)
        U[i,j] = (O[2i-1, 2j-1] + O[2i, 2j] + 1im * (O[2i, 2j-1] - O[2i-1, 2j]))/2
    end
    return U
end

"""
    adjoint(ψ::MajoranaReg{T}) * φ::MajoranaReg{T} -> Complex{T}

Compute the inner product ``⟨ψ|φ⟩`` between two `MajoranaReg` states with definite occupation.

!!! warning
    Both `ψ` and `φ` must be states with a definite occupation (i.e. prepared
    by a number-conserving FLO circuit from a product state). This function will error
    if either register does not have a definite particle number.

!!! warning
    The result is only consistent with the `ArrayReg` behaviour if the number-conserving
    unitary `U` that prepared the state leaves the vacuum / all-zero state invariant,
    i.e. ``U|Ω⟩ = |Ω⟩`` without any relative phases. Otherwise, the result will differ
    from the `ArrayReg` overlap by that phase.
"""
function Base.:*(bra::AdjointMajoranaReg, ket::MajoranaReg)
    nqubits(bra) == nqubits(ket) || error(
        "expect ⟨bra| and |ket⟩ to have the same number of qubits, got $(nqubits(bra)) and $(nqubits(ket))"
    )
    ψ = parent(bra)
    xlocs1, R1 = productstate_and_circuit_decomposition(ψ)
    xlocs2, R2 = productstate_and_circuit_decomposition(ket)
    length(xlocs1) == length(xlocs2) || return zero(complex(eltype(ψ)))
    S = majorana2unitary(R2' * R1)
    n_swaps1 = sum(xlocs1) - length(xlocs1) * (length(xlocs1) - 1) ÷ 2
    n_swaps2 = sum(xlocs2) - length(xlocs2) * (length(xlocs2) - 1) ÷ 2
    permutation_sign = (-1)^(n_swaps1 + n_swaps2)
    return permutation_sign * det(S[xlocs2, xlocs1])
end
