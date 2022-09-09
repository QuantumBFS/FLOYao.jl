# ------------------------------
# Calculating expectation values
# ------------------------------

"""
    majorana_expect(block::AbstractMatrix, locs, reg::MajoranaReg)

Calculate `⟨reg|block|reg⟩` where `block` is a hamiltonian written as the 
coefficients in a sum of squares of Majorana operators and `reg` a 
MajoranaReg.
"""
function majorana_expect(block::AbstractMatrix, locs, reg::MajoranaReg)
    fullH = zero(reg.state)
    matlocs = (2*(locs[1]-1)+1:2(locs[end]))
    fullH[matlocs,matlocs] .= block
    
    expval = sum(1:nqubits(reg), init=zero(eltype(reg.state))) do i
        reg.state[:,2i] ⋅ (fullH * reg.state[:,2i-1]) - reg.state[:,2i-1] ⋅ (fullH * reg.state[:,2i])
    end
    
    offset_expval = sum(1:nqubits(reg), init=zero(eltype(reg.state))) do i
        - reg.state[:,2i] ⋅ (fullH * reg.state[:,2i]) - reg.state[:,2i-1] ⋅ (fullH * reg.state[:,2i-i])
    end |> imag
    
    return (- expval - offset_expval) / 4
end

for BT in [:AbstractAdd, :Scale, :(KronBlock{2}), :(PutBlock{2,1,ZGate})]
    @eval function Yao.expect(op::$BT, reg::MajoranaReg)
        YaoBlocks._check_size(reg, op)
        majoranaham = yaoham2majoranasquares(eltype(reg), op)
        return majorana_expect(majoranaham, 1:nqubits(reg), reg)
    end
end

Yao.expect(op::Scale, reg::MajoranaReg) = op.alpha * Yao.expect(op.content, reg)
#
# ----------------------
# Calculating fidelities
# ----------------------
"""
    bitstring_probability(reg::MajoranaReg, bit_string::BitStr)

The probability to measure a given `bit_string` when measuring `reg`
"""
function bitstring_probability(reg::MajoranaReg{T}, bit_string::DitStr{2,N,ST}) where {T,N,ST}
    @assert nqubits(reg) == N
    p = one(T)
    M = covariance_matrix(reg)
    for i in 1:N
        ni = bit_string[i]
        pi = (1 + (-1)^ni * M[2i-1,2i]) / 2
        pi ≈ 0 && return zero(T)
        p *= pi
        update_covariance_matrix!(M, i, pi, ni)
    end
    return p
end

function Yao.fidelity(reg1::MajoranaReg{T1}, reg2::MajoranaReg{T2}) where {T1,T2}
    T = promote_type(T1, T2)
    nq = nqubits(reg1)
    reg = FLOYao.zero_state(T, nq)
    reg.state = reg2.state' * reg1.state
    return √bitstring_probability(reg, BitStr{nq}(0))
end
