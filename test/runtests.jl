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

using Test
using FLOYao
using Yao.YaoAPI # needed for the doctest tests
using Yao
using StatsBase
using Documenter
using LinearAlgebra
import FLOYao: majorana2arrayreg, NonFLOException

# Run the doctests
DocMeta.setdocmeta!(FLOYao, :DocTestSetup, :(using FLOYao, Yao); recursive=true)
doctest(FLOYao)

@const_gate TestGate::ComplexF64 = [1 0 ; 0 exp(-1im*π/5)]
@const_gate FSWAP::ComplexF64 = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 -1]
@const_gate ManualX::ComplexF64 = [0 1; 1 0]

function ising_hamiltonian(nq, J, h)
    U = map(1:nq-1) do i
        J * kron(nq, i => X, i+1 => X)
    end |> sum

    T = map(1:nq) do i
        h * kron(nq, i => Z)
    end |> sum

    hamiltonian = T + U

    return hamiltonian
end

@testset "fast_overlap" begin
    nq = 4
    x = randn(ComplexF64, 10, 10)
    y = randn(ComplexF64, 10, 10)
    dmat = [rand() < 0.4 ? randn() : 0.0 for i=1:10, j=1:10]
    B = FLOYao.SparseMatrixCOO(dmat)
    r2 = FLOYao.fast_overlap(y, B, x)
    @test isapprox(r2, tr(y' * B * x), atol=1e-7)
    @test isapprox(r2, y ⋅ (B * x), atol=1e-7)
end

@testset "fast_add!" begin
    A = randn(10, 10)
    dmat = [rand() < 0.4 ? randn() : 0.0 for i=1:10, j=1:10]
    D = FLOYao.SparseMatrixCOO(dmat)

    C = copy(A)
    @test FLOYao.fast_add!(C, D) ≈ A + D
    C = copy(A)
    @test FLOYao.fast_add!(C, Matrix(D)) ≈ A + D

end

@testset "MajoranaRegister" begin
    nq = 2
    mreg = FLOYao.zero_state(nq)
    areg = zero_state(nq)
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    FLOYao.one_state!(mreg)
    areg = product_state(bit"11")
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    mreg = FLOYao.product_state(bit"10101")
    areg = product_state(bit"10101")
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    mreg2 = FLOYao.product_state(bit"11111")
    FLOYao.one_state!(mreg)
    @test fidelity(mreg, mreg2) ≈ 1.
    @test fidelity(mreg, FLOYao.one_state(nqubits(mreg))) ≈ 1.

    r1 = FLOYao.product_state(Float32, bit"001")
    r2 = FLOYao.product_state(Float32, [1, 0, 0])
    r3 = FLOYao.product_state(Float32, 3, 0b001)
    @test r1 ≈ r2   # because we read bit strings from right to left, vectors from left to right.
    @test r1 ≈ r3

    nq = 5
    mreg = FLOYao.rand_state(nq)
    @test mreg.state * mreg.state' ≈ I(2nq)
    mreg = FLOYao.rand_state(Float32, nq)
    @test mreg.state * mreg.state' ≈ I(2nq)

    tmp = copy(mreg)
    @test state(tmp) ≈ state(tmp)

    tmp = similar(mreg)
    copyto!(tmp, mreg)
    @test state(tmp) ≈ state(tmp)

    tmp_wrongsize = FLOYao.zero_state(nqubits(mreg) + 1)
    @test_throws DimensionMismatch copyto!(tmp_wrongsize, mreg)
end

@testset "PutBlock" begin
    nq = 3
    θ = 0.4
    mreg = FLOYao.zero_state(nq)
    areg = zero_state(nq)
    xxg = kron(nq, 1 => X, 2 => X)

    rg = rot(xxg, θ)
    mreg |> rg
    areg |> rg
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    pb = put(nq, 2 => Rz(θ))
    mreg |> pb
    areg |> pb
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    xg = put(nq, 2 => X)
    mreg |> xg
    areg |> xg
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    rz = put(nq, 2 => Rz(θ))
    mreg |> rz
    areg |> rz
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    fswap = put(nq, (2, 3) => FSWAP)
    @test_warn "Calling manual instruct!" mreg |> fswap
    areg |> fswap
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    yg = put(nq, 3 => Y)
    mreg |> yg
    areg |> yg
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    tg = put(nq, 2 => T)
    mreg |> tg
    areg |> tg
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    tdg = put(nq, 2 => T')
    mreg |> tdg
    areg |> tdg
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    @test_warn "Calling manual instruct!" Yao.instruct!(mreg, mat(FSWAP), (1,2))
    Yao.instruct!(areg, mat(FSWAP), (1,2))
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    rx = put(nq, 2 => Rx(0.1))
    ry = put(nq, 1 => Ry(0.1))
    @test_throws NonFLOException mreg |> rx
    @test_throws NonFLOException mreg |> ry
end

@testset "KronBlock" begin
    nq = 3
    θ = 0.4
    mreg = FLOYao.zero_state(nq)
    areg = zero_state(nq)
    xxg = kron(nq, 1 => X, 2 => X)
    rg = rot(xxg, θ)

    mreg |> rg
    areg |> rg
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    mreg |> xxg
    areg |> xxg
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    kb = kron(nq, 2 => Rz(0.4), 3 => shift(1.))
    mreg |> kb
    areg |> kb
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    xx = kron(nq, 1 => ManualX, 2 => ManualX)
    @test_warn "Calling manual instruct!" mreg |> xx
    areg |> xx
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    xyrot = kron(nq, 2:3 => rot(kron(2, 1=>X, 2=>Y), θ))
    mreg |> xyrot
    areg |> xyrot
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.
end

@testset "RepeatedBlock" begin
    nq = 3
    θ = 0.4
    mreg = FLOYao.zero_state(nq)
    areg = zero_state(nq)
    xxg = kron(nq, 1 => X, 2 => X)
    rg = rot(xxg, θ)
    mreg |> rg
    areg |> rg

    rp = repeat(nq, X, 1:2)
    mreg |> rp
    areg |> rp
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    rp = repeat(nq, TestGate, 1:2)
    @test_warn "Calling manual instruct!" mreg |> rp
    areg |> rp
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.
end

@testset "TimeEvolution" begin
    nq = 5
    J = 1.5
    h = -1.
    hamiltonian = ising_hamiltonian(nq, J, h)
    ising_evolution = time_evolve(hamiltonian, 1.)

    mreg = FLOYao.zero_state(nq)
    areg = zero_state(nq)

    mreg |> ising_evolution
    areg |> ising_evolution
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.
end

@testset "expect" begin
    nq = 4
    circuit = chain(nq)

    θ = π/8
    xxg = kron(nq, 1 => X, 2 => Y)
    rg = rot(xxg, θ)
    push!(circuit, rg)
    push!(circuit, put(nq, 3=>Rz(0.5)))
    push!(circuit, rg)

    θ = π/5
    xxg = kron(nq, 2 => X, 3 => Z, 4 => Y)
    rg = rot(xxg, θ)
    rz = put(nq, 3 => Rz(θ))
    push!(circuit, rg)
    push!(circuit, put(nq, 3=>Rz(0.5)))
    push!(circuit, put(nq, (2,3) => FSWAP))
    push!(circuit, put(nq, 1=>Z))
    push!(circuit, put(nq, 4=>X))
    push!(circuit, rg)
    push!(circuit, rz)

    ham = (put(nq, 1=>Z) + 2kron(nq, 1=>X, 2=>Z, 3=>Z, 4=>X) + 3.5put(nq, 2=>Z)
          + 0.5kron(nq, 1=>Z, 2=>Z) - kron(nq, 2 => X, 4 => Y) # some higher order terms
          + 1.5kron(nq, 1=>X) - kron(nq, 1=>Z, 3=>Y)           # some odd parity terms
          )
    mreg = FLOYao.zero_state(nq)
    areg = zero_state(nq)
    mreg |> put(nq, 2=>X)
    areg |> put(nq, 2=>X)

    meval = expect(ham, mreg |> circuit)
    aeval = expect(ham, areg |> circuit)
    @test meval ≈ aeval

    meval = expect(ham[end], mreg)
    aeval = expect(ham[end], areg)
    @test meval ≈ aeval atol=1e-15

    op = kron(nq, 2=>X, 3=>X)
    meval = expect(op, mreg)
    aeval = expect(op, areg)
    @test meval ≈ aeval

    mreg = FLOYao.rand_state(nq)
    areg = FLOYao.majorana2arrayreg(mreg)
    ops = [
        kron(nq, 1=>X, 2=>X)
        kron(nq, 2=>X, 3=>X)
        kron(nq, 1=>X, 3=>X)
        kron(nq, 1=>X, 4=>X)
        kron(nq, 1=>Y, 2=>Y)
        kron(nq, 2=>Y, 3=>Y)
        kron(nq, 1=>Y, 3=>Y)
        kron(nq, 1=>Y, 4=>Y)
        kron(nq, 1=>X, 2=>Y)
        kron(nq, 2=>X, 3=>Y)
        kron(nq, 1=>X, 3=>Y)
        kron(nq, 1=>X, 4=>Y)
        kron(nq, 1=>Y, 2=>X)
        kron(nq, 2=>Y, 3=>X)
        kron(nq, 1=>Y, 3=>X)
        kron(nq, 1=>Y, 4=>X)
        kron(nq, 1=>X, 2=>Z, 3=>X)
        kron(nq, 1=>X, 2=>Z, 3=>Z, 4=>X)
        kron(nq, 1=>X, 2=>Z, 3=>Y)
        kron(nq, 1=>X, 2=>Z, 3=>Z, 4=>Y)
        kron(nq, 1=>Z)
        kron(nq, 2=>Z)
        kron(nq, 2=>Z, 3=>Z)
        kron(nq, 1=>Z, 3=>Z)
        kron(nq, 1=>Z, 4=>Z)
    ]
    for op in ops
        meval = expect(op, mreg)
        aeval = expect(op, areg)
        @test meval ≈ aeval
    end

    # test kronblocks inside put blocks
    ham = put(nq, 2:3 => kron(2, 1 => X, 2 => Y)) + 2put(nq, 1:2 => kron(2, 1 => X, 2 => X))
    meval = expect(ham, mreg)
    aeval = expect(ham, areg)
    @test meval ≈ aeval
end

@testset "autodiff" begin
    nq = 4
    circuit = chain(nq)

    θ = π/8
    xxg = kron(nq, 1 => X, 2 => Z, 3 => X)
    rg = rot(xxg, θ)
    push!(circuit, rg)
    push!(circuit, put(nq, 2=>X))
    push!(circuit, put(nq, 3=>Rz(0.5)))
    push!(circuit, put(nq, (2,3) => FSWAP))
    push!(circuit, rg)

    θ = π/5
    xxg = kron(nq, 2 => X, 3 => Z, 4 => Y)
    rg = rot(xxg, θ)
    push!(circuit, rg)
    push!(circuit, put(nq, 3=>Rz(0.5)))
    push!(circuit, put(nq, 1=>Z))

    J = 1.5
    h = -1.
    hamiltonian = ising_hamiltonian(nq, J, h)
    ising_evolution = time_evolve(hamiltonian, 1.)
    push!(circuit, ising_evolution)

    ham = (put(nq, 1=>Z) + 2kron(nq, 1=>X, 2=>Z, 3=>Z, 4=>X) + 3.5put(nq, 2=>Z)
          + 0.5kron(nq, 1=>Z, 2=>Z) - kron(nq, 2 => X, 4 => Y)
          )
    majoranasum = MajoranaSum(ham)
    mreg = FLOYao.zero_state(nq)
    areg = zero_state(nq)
    mreg |> put(nq, 1=>X)
    areg |> put(nq, 1=>X)
    mgrad1 = @test_warn "Calling manual" expect'(ham, mreg => circuit)[2]
    mgrad2 = @test_warn "Calling manual" expect'(majoranasum, mreg => circuit)[2]
    agrad = expect'(ham, areg => circuit)[2]

    @test mgrad1 ≈ agrad
    @test mgrad2 ≈ agrad

    mreg = FLOYao.zero_state(nq)
    areg = zero_state(nq)
    mreg |> put(nq, 1=>X) |> circuit
    areg |> put(nq, 1=>X) |> circuit

    mreg1 = copy(mreg)
    mreg2 = copy(mreg)

    mregδ1 = Yao.AD.expect_g(ham, mreg1)
    mregδ2 = Yao.AD.expect_g(majoranasum, mreg2)
    aregδ = Yao.AD.expect_g(ham, areg)

    (in_mreg1, in_mregδ1), params_mregδ1 = @test_warn "Calling manual" Yao.AD.apply_back((mreg1, mregδ1), circuit)
    (in_mreg2, in_mregδ2), params_mregδ2 = @test_warn "Calling manual" Yao.AD.apply_back((mreg2, mregδ2), circuit)
    (in_areg, in_aregδ), params_aregδ = Yao.AD.apply_back((areg, aregδ), circuit)
    @test params_mregδ1 ≈ params_aregδ
    @test params_mregδ2 ≈ params_aregδ

    @test fidelity(majorana2arrayreg(in_mreg1), in_areg) ≈ 1.
    @test fidelity(majorana2arrayreg(in_mreg2), in_areg) ≈ 1.

    mzero = FLOYao.zero_state(nq)
    azero = zero_state(nq)
    mreg = apply(mzero, put(nq, 1=>X))
    areg = apply(azero, put(nq, 1=>X))
    mgrad1 = @test_warn "Calling manual" fidelity'(mzero, mreg => circuit)[2][2]
    agrad = fidelity'(azero, areg => circuit)[2][2]
    @test mgrad1 ≈ agrad

    mreg2 = FLOYao.rand_state(nq)
    if det(mreg2.state) < 0
        apply!(mreg2, put(nq, 1=>X))
    end
    areg2 = majorana2arrayreg(mreg2)
    mgrad1 = @test_warn "Calling manual" fidelity'(mreg => circuit, mreg2)[1][2]
    agrad = fidelity'(areg => circuit, areg2)[1][2]
    @test mgrad1 ≈ agrad

    # testing singular matrices
    mreg = FLOYao.zero_state(nq)
    areg = zero_state(nq)

    xxg = kron(nq, 1 => X, 2 => X)
    circuit = rot(xxg, π/2)

    ham = kron(nq, 1=>Z)
    agrad = expect'(ham, areg => circuit)[2]
    mgrad = expect'(ham, mreg => circuit)[2]
    @test agrad ≈ mgrad

    circuit = rot(xxg, 0.)

    ham = kron(nq, 1=>X, 2=>X)
    agrad = expect'(ham, areg => circuit)[2]
    mgrad = expect'(ham, mreg => circuit)[2]
    @test agrad ≈ mgrad
end

@testset "fidelity" begin
    nq = 5
    mreg1 = FLOYao.rand_state(nq)
    areg1 = FLOYao.majorana2arrayreg(mreg1)
    mreg2 = FLOYao.rand_state(nq)
    areg2 = FLOYao.majorana2arrayreg(mreg2)

    @test isapprox(fidelity(mreg1, mreg2), fidelity(areg1, areg2), atol=1e-7)
    @test isapprox(fidelity(mreg1, FLOYao.zero_state_like(mreg1)),
                   fidelity(areg1, zero_state_like(areg1, nq)),
                   atol=1e-7)

    # ensure that if previously the fidelity was zero because of different
    # determinants it now isn't
    x1 = put(nq, 1 => X)
    mreg1 |> x1
    areg1 |> x1
    @test isapprox(fidelity(mreg1, mreg2), fidelity(areg1, areg2), atol=1e-7)
    @test isapprox(fidelity(mreg1, FLOYao.zero_state_like(mreg1)),
                   fidelity(areg1, zero_state_like(areg1, nq)),
                   atol=1e-7)
end

@testset "measure" begin
    nq = 4
    circuit = chain(nq)

    θ = π/8
    xxg = kron(nq, 1 => X, 2 => Y)
    rg = rot(xxg, θ)
    push!(circuit, rg)
    push!(circuit, put(nq, 3=>Rz(0.5)))
    push!(circuit, rg)

    θ = π/5
    xxg = kron(nq, 2 => X, 3 => Z, 4 => Y)
    rg = rot(xxg, θ)
    push!(circuit, rg)
    push!(circuit, put(nq, 3=>Rz(0.5)))
    push!(circuit, put(nq, 1=>Z))
    push!(circuit, put(nq, 4=>X))
    push!(circuit, rg)

    mreg = FLOYao.product_state(bit"1001")
    areg = product_state(bit"1001")
    mreg |> circuit
    areg |> circuit

    testreg = copy(mreg)
    res = measure!(testreg)
    @test fidelity(testreg, FLOYao.product_state(res)) ≈ 1.

    testreg = copy(mreg)
    res = measure!(NoPostProcess(), testreg)
    @test fidelity(testreg, FLOYao.product_state(res)) ≈ 1.

    testreg = copy(mreg)
    bits = bit"1100"
    res = measure!(ResetTo(bits), testreg)
    @test fidelity(testreg, FLOYao.product_state(bits)) ≈ 1.

    mprobs = [FLOYao.bitstring_probability(mreg, BitStr{nq}(b)) for b in 0:2^nq-1]
    aprobs = abs2.(areg.state)
    @test mprobs ≈ aprobs

    msamples = measure(mreg, nshots=100000)
    asamples = measure(areg, nshots=100000)
    mhist = StatsBase.normalize(fit(Histogram, Int.(msamples), nbins=2^nq), mode=:probability)
    ahist = StatsBase.normalize(fit(Histogram, Int.(asamples), nbins=2^nq), mode=:probability)
    @test sum(abs, ahist.weights - mhist.weights) < 0.01


    msamples = measure(mreg, [3,1,4], nshots=100000)
    asamples = measure(areg, [3,1,4], nshots=100000)
    mhist = StatsBase.normalize(fit(Histogram, Int.(msamples), nbins=2^3), mode=:probability)
    ahist = StatsBase.normalize(fit(Histogram, Int.(asamples), nbins=2^3), mode=:probability)
    @test sum(abs, ahist.weights - mhist.weights) < 0.01
end

@testset "utils" begin
    nq = 4
    ham = put(nq, 1=>Z) + 2kron(nq, 1=>X, 2=>Z, 3=>Z, 4=>X) + 3.5put(nq, 2=>Z)
    ham_qubit = mat(ham)
    ham_pauli = FLOYao.qubit2paulibasis(ham_qubit)
    @test FLOYao.paulibasis2qubitop(ham_pauli) ≈ ham_qubit

    ham_majorana = FLOYao.paulibasis2majoranasquares(ham_pauli)
    @test FLOYao.majoranasquares2qubitbasis(ham_majorana) ≈ ham_qubit
end

@testset "docs" begin
    doctest(FLOYao; manual=false)
end
