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

using Pkg
Pkg.develop(path="../../..")    # Add FLOYao as development dependency
Pkg.develop("YaoArrayRegister") # Needed until 0.9.10 is released

using Test
using FLOYao
using Yao
using CUDA
using StatsBase
using LinearAlgebra
using Random

import FLOYao: majorana2arrayreg, NonFLOException
import ExponentialUtilities    # needed to trigger loading of FLOYaoCUDAExt

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

@testset "CuMajoranaRegister" begin
    nq = 2
    mreg = FLOYao.cuzero_state(nq)
    areg = cuzero_state(ComplexF32, nq)
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ one(eltype(mreg))

    mreg2 = FLOYao.zero_state(nq) |> cu
    @test mreg2 == mreg

    FLOYao.one_state!(mreg)
    areg = cuproduct_state(bit"11")
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ one(eltype(mreg))

    mreg = FLOYao.cuproduct_state(bit"10101")
    areg = cuproduct_state(bit"10101")
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ one(eltype(mreg))

    mreg2 = FLOYao.cuproduct_state(bit"11111")
    FLOYao.one_state!(mreg)
    @test fidelity(mreg, mreg2) ≈ one(eltype(mreg))
    @test fidelity(mreg, FLOYao.cuone_state(nqubits(mreg))) ≈ one(eltype(mreg))

    r1 = FLOYao.cuproduct_state(Float32, bit"001")
    r2 = FLOYao.cuproduct_state(Float32, [1, 0, 0])
    r3 = FLOYao.cuproduct_state(Float32, 3, 0b001)
    @test r1 ≈ r2   # because we read bit strings from right to left, vectors from left to right.
    @test r1 ≈ r3

    nq = 5
    mreg = FLOYao.rand_state(nq) |> cu
    @test mreg.state * mreg.state' ≈ CUDA.CuMatrix(I(2nq))

    tmp = copy(mreg)
    @test state(tmp) ≈ state(tmp)

    tmp = similar(mreg)
    copyto!(tmp, mreg)
    @test state(tmp) ≈ state(tmp)
end

@testset "PutBlock" begin
    nq = 3
    θ = 0.4
    mreg = FLOYao.cuzero_state(nq)
    areg = cuzero_state(nq)
    xxg = kron(nq, 1 => X, 2 => X)

    rg = rot(xxg, θ)
    mreg |> rg
    areg |> rg
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ one(eltype(mreg))

    pb = put(nq, 2 => Rz(θ))
    mreg |> pb
    areg |> pb
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ one(eltype(mreg))

    xg = put(nq, 2 => X)
    mreg |> xg
    areg |> xg
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ one(eltype(mreg))

    rz = put(nq, 2 => Rz(θ))
    mreg |> rz
    areg |> rz
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ one(eltype(mreg))

    yg = put(nq, 3 => Y)
    mreg |> yg
    areg |> yg
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ one(eltype(mreg))

    tg = put(nq, 2 => T)
    mreg |> tg
    areg |> tg
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ one(eltype(mreg))

    tdg = put(nq, 2 => T')
    mreg |> tdg
    areg |> tdg
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ one(eltype(mreg))
end

@testset "KronBlock" begin
    nq = 3
    θ = 0.4
    mreg = FLOYao.cuzero_state(nq)
    areg = cuzero_state(nq)
    xxg = kron(nq, 1 => X, 2 => X)
    rg = rot(xxg, θ)

    mreg |> rg
    areg |> rg
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ one(eltype(mreg))

    mreg |> xxg
    areg |> xxg
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ one(eltype(mreg))

    kb = kron(nq, 2 => Rz(0.4), 3 => shift(1.))
    mreg |> kb
    areg |> kb
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ one(eltype(mreg))

    xyrot = kron(nq, 2:3 => rot(kron(2, 1=>X, 2=>Y), θ))
    mreg |> xyrot
    areg |> xyrot
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ one(eltype(mreg))
end

@testset "RepeatedBlock" begin
    nq = 3
    θ = 0.4
    mreg = FLOYao.cuzero_state(nq)
    areg = cuzero_state(nq)
    xxg = kron(nq, 1 => X, 2 => X)
    rg = rot(xxg, θ)
    mreg |> rg
    areg |> rg

    rp = repeat(nq, X, 1:2)
    mreg |> rp
    areg |> rp
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ one(eltype(mreg))
end

@testset "TimeEvolution" begin
    nq = 5 
    J = 1.5 
    h = -1.
    hamiltonian = ising_hamiltonian(nq, J, h)
    ising_evolution = time_evolve(hamiltonian, 1.)

    mreg = FLOYao.cuzero_state(nq)
    areg = cuzero_state(nq)

    mreg |> ising_evolution
    areg |> ising_evolution
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ one(eltype(mreg))
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
    push!(circuit, put(nq, 1=>Z))
    push!(circuit, put(nq, 4=>X))
    push!(circuit, rg)  
    push!(circuit, rz)  

    ham = put(nq, 1=>Z) + 2kron(nq, 1=>X, 2=>Z, 3=>Z, 4=>X) + 3.5f0put(nq, 2=>Z)
    mreg = FLOYao.cuzero_state(nq)
    areg = cuzero_state(nq)
    mreg |> put(nq, 2=>X)
    areg |> put(nq, 2=>X)

    meval = expect(ham, mreg |> circuit)
    aeval = expect(ham, areg |> circuit) 
    @test meval ≈ aeval

    meval = expect(ham[end], mreg)
    aeval = expect(ham[end], areg)
    @test meval ≈ aeval

    # test kronblocks inside put blocks
    ham = put(nq, 2:3 => kron(2, 1 => X, 2 => Y)) + 2put(nq, 1:2 => kron(2, 1 => X, 2 => Y))
    meval = expect(ham, mreg)
    aeval = expect(ham, areg)
    @test meval ≈ aeval atol=1e-7
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

    ham = put(nq, 1=>Z) + 2kron(nq, 1=>X, 2=>Z, 3=>Z, 4=>X) + 3.5put(nq, 2=>Z)
    mreg = FLOYao.cuzero_state(nq)
    areg = cuzero_state(nq)
    mreg |> put(nq, 1=>X)
    areg |> put(nq, 1=>X)
    mgrad = expect'(ham, mreg => circuit)[2]
    agrad = expect'(ham, areg => circuit)[2]

    @test mgrad ≈ agrad

    mreg = FLOYao.cuzero_state(nq)
    areg = cuzero_state(nq)
    mreg |> put(nq, 1=>X) |> circuit
    areg |> put(nq, 1=>X) |> circuit

    mregδ = Yao.AD.expect_g(ham, mreg)
    aregδ = Yao.AD.expect_g(ham, areg)

    (in_mreg, in_mregδ), params_mregδ = Yao.AD.apply_back((mreg, mregδ), circuit)
    (in_areg, in_aregδ), params_aregδ = Yao.AD.apply_back((areg, aregδ), circuit)
    @test params_mregδ ≈ params_aregδ

    @test fidelity(majorana2arrayreg(in_mreg), in_areg) ≈ one(eltype(in_mreg))
end

@testset "fidelity" begin
    nq = 5
    Random.seed!(42)
    mreg1 = FLOYao.rand_state(nq) |> cu
    areg1 = FLOYao.majorana2arrayreg(mreg1)
    mreg2 = FLOYao.rand_state(nq) |> cu
    areg2 = FLOYao.majorana2arrayreg(mreg2)

    @test isapprox(fidelity(mreg1, mreg2), fidelity(areg1, areg2), atol=1e-5)
    @test isapprox(fidelity(mreg1, FLOYao.zero_state_like(mreg1) |> cu),
                   fidelity(areg1, zero_state_like(areg1, nq) |> cu),
                   atol=1e-5)

    # ensure that if previously the fidelity was zero because of different 
    # determinants it now isn't
    x1 = put(nq, 1 => X)
    mreg1 |> x1
    areg1 |> x1
    @test isapprox(fidelity(mreg1, mreg2), fidelity(areg1, areg2), atol=1e-5)
    @test isapprox(fidelity(mreg1, FLOYao.zero_state_like(mreg1) |> cu),
                   fidelity(areg1, zero_state_like(areg1, nq)),
                   atol=1e-5)
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

    mreg = FLOYao.cuproduct_state(bit"1001")
    areg = cuproduct_state(bit"1001")
    mreg |> circuit
    areg |> circuit

    testreg = copy(mreg)
    res = measure!(testreg)
    @test fidelity(testreg, FLOYao.cuproduct_state(res)) ≈ one(eltype(testreg))

    testreg = copy(mreg)
    res = measure!(NoPostProcess(), testreg)
    @test fidelity(testreg, FLOYao.cuproduct_state(res)) ≈ one(eltype(testreg))

    testreg = copy(mreg)
    bits = bit"1100"
    res = measure!(ResetTo(bits), testreg)
    @test fidelity(testreg, FLOYao.cuproduct_state(bits)) ≈ one(eltype(testreg))

    mprobs = [FLOYao.bitstring_probability(mreg, BitStr{nq}(b)) for b in 0:2^nq-1]
    aprobs = abs2.(areg.state) |> Array
    @test mprobs ≈ aprobs

    msamples = measure(mreg, nshots=100000)
    asamples = measure(areg, nshots=100000)
    mhist = StatsBase.normalize(fit(Histogram, Int.(msamples), nbins=2^nq), mode=:probability)
    ahist = StatsBase.normalize(fit(Histogram, Int.(asamples), nbins=2^nq), mode=:probability)
    @test norm(ahist.weights - mhist.weights, 1) < 0.01

    msamples = measure(mreg, [3,1,4], nshots=100000)
    asamples = measure(areg, [3,1,4], nshots=100000)
    mhist = StatsBase.normalize(fit(Histogram, Int.(msamples), nbins=2^3), mode=:probability)
    ahist = StatsBase.normalize(fit(Histogram, Int.(asamples), nbins=2^3), mode=:probability)
    @test norm(ahist.weights - mhist.weights, 1) < 0.01
end
