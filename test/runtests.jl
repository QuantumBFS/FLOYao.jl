using Yao
using Test
using FLOYao
import FLOYao: majorana2arrayreg

@const_gate TestGate::ComplexF64 = [1 0 ; 0 exp(-1im*π/5)]
@const_gate FSWAP::ComplexF64 = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 -1]

@testset "MajoranaRegister" begin
    nq = 2
    mreg = FLOYao.zero_state(nq)
    areg = zero_state(nq)
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    mreg = FLOYao.zero_state(nq)
    areg = zero_state(nq)
    areg |> kron(2, 1=>X, 2=>X)
    mreg |> kron(2, 1=>X, 2=>X)
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.
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
    push!(circuit, rg)  
    push!(circuit, put(nq, 3=>Rz(0.5)))
    push!(circuit, put(nq, (2,3) => FSWAP))
    push!(circuit, put(nq, 1=>Z))
    push!(circuit, put(nq, 4=>X))
    push!(circuit, rg)  

    ham = put(nq, 1=>Z) + 2kron(nq, 1=>X, 2=>Z, 3=>Z, 4=>X) + 3.5put(nq, 2=>Z)
    mreg = FLOYao.zero_state(nq)
    areg = zero_state(nq)
    mreg |> put(nq, 2=>X)
    areg |> put(nq, 2=>X)

    meval = @test_warn "Calling manual" expect(ham, mreg => circuit)
    aeval = expect(ham, areg => circuit) 

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

    ham = put(nq, 1=>Z) + 2kron(nq, 1=>X, 2=>Z, 3=>Z, 4=>X) + 3.5put(nq, 2=>Z)
    mreg = FLOYao.zero_state(nq)
    areg = zero_state(nq)
    mreg |> put(nq, 1=>X)
    areg |> put(nq, 1=>X)
    mgrad = @test_warn "Calling manual" expect'(ham, mreg => circuit)[2] 
    agrad = expect'(ham, areg => circuit)[2] 

    @test mgrad ≈ agrad
end
