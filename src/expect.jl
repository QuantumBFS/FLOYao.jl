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

for BT in [:AbstractAdd, :(KronBlock{2}), :(PutBlock{2,1,ZGate})]
    @eval function Yao.expect(op::$BT, reg::MajoranaReg)
        YaoBlocks._check_size(reg, op)
        majoranaham = yaoham2majoranasquares(eltype(reg), op)
        return majorana_expect(majoranaham, 1:nqubits(reg), reg)
    end
end

Yao.expect(op::Scale, reg::MajoranaReg) = op.alpha * Yao.expect(op.content, reg)
