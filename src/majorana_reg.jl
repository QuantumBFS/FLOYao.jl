# ------------------------------
# The actual Register Definition
# ------------------------------

mutable struct MajoranaReg{T<:Real} <: AbstractRegister{2}
    state::Matrix{T}
end

"""
    MajoranaReg([T=Float64,] n)

Create a Majorana register on `n` qubits in the vacuum state ``|Ω⟩`` with
storage type `T`.
"""
function MajoranaReg(::Type{T}, n::Integer) where {T}
    return MajoranaReg(Matrix{T}(I(2n)))
end

MajoranaReg(n::Integer) = MajoranaReg(Float64, n)
Yao.nqubits(reg::MajoranaReg) = size(reg.state, 1) ÷ 2
Yao.nqudits(reg::MajoranaReg) = size(reg.state, 1) ÷ 2
Yao.nactive(reg::MajoranaReg) = Yao.nqubits(reg)
Yao.nbatch(reg::MajoranaReg) = 1
Yao.nremain(reg::MajoranaReg) = 0
Base.copy(reg::MajoranaReg) = MajoranaReg(copy(reg.state))
Base.eltype(reg::MajoranaReg) = eltype(reg.state)

# The detailed version showing the contents of the register in e.g. 
# the jupyter cell output
function Base.show(io::IO, m::MIME"text/plain", reg::MajoranaReg)
    println(io, typeof(reg), " with $(Yao.nqubits(reg)) qubits:")
    show(io, m, reg.state)
end

# Less detailed version that is used e.g. in string interpolations
function Base.show(io::IO, reg::MajoranaReg)
    print(io, typeof(reg), "($(nqubits(reg)))")
end

"""
    majorana2arrayreg(reg::MajoranaReg)

Converts a `2n×2n` MajoranaReg `reg` into a `2^n` ArrayReg.

# Note
This implementation is not very clever and should mainly be used for debugging
purposes with small numbers of qubits. If ⟨Ω|U|Ω⟩ is close to zero, it is not
very accurate.
"""
function majorana2arrayreg(reg::MajoranaReg)
    nq = nqubits(reg)

    # praying here, that the zero_state has non-zero overlap with the 
    # new vacuum state.
    # This first bit gets areg into the new vacuum by piping it through
    # the projector U|Ω⟩⟨Ω|U^† 
    areg = uniform_state(Complex{eltype(reg)}, nq)
    for i in 1:nq
        γ_i1 = majoranaop(nq, reg.state[:,2i-1])
        γ_i2 = majoranaop(nq, reg.state[:,2i])
        circuit = 1im*chain(nq, γ_i2, γ_i1) + igate(nq)
        areg |> circuit
    end
    normalize!(areg)
    return areg
end
