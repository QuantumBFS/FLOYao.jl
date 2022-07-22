using Yao
using Test
using FLOYao
import FLOYao: majorana2arrayreg

@testset "MajoranaRegister" begin
    nq = 2
    mreg = MajoranaReg(nq)
    areg = zero_state(nq)
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    mreg = MajoranaReg(nq)
    areg = zero_state(nq)
    areg |> kron(2, 1=>X, 2=>X)
    mreg |> kron(2, 1=>X, 2=>X)
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.
end

@testset "instruct" begin
    nq = 3
    θ = 0.4
    mreg = MajoranaReg(nq)
    areg = zero_state(nq)
    pb = put(nq, 2 => Rz(θ))
    xxg = kron(nq, 1 => X, 2 => X)
    rg = rot(xxg, θ)

    mreg |> pb
    areg |> pb
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    mreg |> put(nq, 2 => X)
    areg |> put(nq, 2 => X)
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    mreg |> rg
    areg |> rg
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.

    mreg |> put(nq, 1 => Y)
    areg |> put(nq, 1 => Y)
    @test fidelity(majorana2arrayreg(mreg), areg) ≈ 1.


    mreg |> xxg
    areg |> xxg
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
    push!(circuit, put(nq, 1=>Z))
    push!(circuit, put(nq, 4=>X))
    push!(circuit, rg)  

    ham = put(nq, 1=>Z) + 2kron(nq, 1=>X, 2=>Z, 3=>Z, 4=>X) + 3.5put(nq, 2=>Z)
    mreg = MajoranaReg(nq)
    areg = zero_state(nq)
    mreg |> put(nq, 2=>X)
    areg |> put(nq, 2=>X)

    @test expect(ham, mreg => circuit) ≈ expect(ham, areg => circuit) 
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

    ham = put(nq, 1=>Z) + 2kron(nq, 1=>X, 2=>Z, 3=>Z, 4=>X) + 3.5put(nq, 2=>Z)
    mreg = MajoranaReg(nq)
    areg = zero_state(nq)
    mreg |> put(nq, 1=>X)
    areg |> put(nq, 1=>X)
    mgrad = expect'(ham, mreg => circuit)[2] 
    agrad = expect'(ham, areg => circuit)[2] 

    @test mgrad ≈ agrad
end
