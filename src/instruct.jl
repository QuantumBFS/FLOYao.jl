# ------------------------------
# Applying gates to the register
# ------------------------------

function Yao.instruct!(reg::MajoranaReg, ::Val{:X}, locs::Tuple)
    loc = locs[1]
    reg.state[2loc:end,:] .*= -1
    return reg
end

function Yao.instruct!(reg::MajoranaReg, ::Val{:Y}, locs::Tuple)
    loc = locs[1]
    reg.state[2loc+1:end,:] .*= -1
    reg.state[2loc-1,:] .*= -1
    return reg
end

function Yao.instruct!(reg::MajoranaReg, ::Val{:Z}, locs::Tuple)
    loc = locs[1]
    reg.state[2loc-1:2loc,:] .*= -1
    return reg
end

for (G, SC) in zip([:T, :Tdag], [(1/√2, 1/√2), (-1/√2, 1/√2)])
    @eval function Yao.instruct!(reg::MajoranaReg, ::Val{$(QuoteNode(G))}, locs::Tuple)
        loc = locs[1]
        s, c = $SC
        ψ1, ψ2 = reg.state[2loc-1,:], reg.state[2loc,:]
        reg.state[2loc-1,:] .= c .* ψ1 .+ s .* ψ2
        reg.state[2loc,:] .= c .* ψ2 .- s .* ψ1
        return reg
    end
end

for G in [:Rz, :Shift]
    @eval function Yao.instruct!(reg::MajoranaReg, ::Val{$(QuoteNode(G))}, locs::Tuple, theta)
        loc = locs[1]
        s, c = sincos(theta)
        ψ1, ψ2 = reg.state[2loc-1,:], reg.state[2loc,:]
        reg.state[2loc-1,:] .= c .* ψ1 .+ s .* ψ2
        reg.state[2loc,:] .= c .* ψ2 .- s .* ψ1
        return reg
    end
end

# Fallback for generic matrix gates. This is _much_ slower than the explicitely
# defined previous gates, so you do profit from defining unsafe_apply! and 
# instruct! for your fancy new FLO gate like I did before.
function Yao.instruct!(reg::MajoranaReg, gate::AbstractMatrix, locs)
    # unless the gate is a tensor product of Z's (but _not_ using Yaos kron
    # block a gate acting on non-consecutive qubits will never be a 
    # FLO gate
    areconsecutive(locs) || throw(NonFLOException("$gate on $locs is not a FLO gate"))
    W = qubit2majoranaevolution(gate, locs)
    matlocs = (2*(locs[1]-1)+1:2(locs[end]))
    @warn """Calling manual instruct!($reg, $gate, $locs).
    You can greatly speed up your FLO gates by exactly implementing unsafe_apply!()
    and instruct!() for them. See FLOYao/src/instruct.jl and  FLOYao/src/apply_composite.jl
    for how to do that.
    """
    reg.state[matlocs,:] .= W * reg.state[matlocs,:]
end

