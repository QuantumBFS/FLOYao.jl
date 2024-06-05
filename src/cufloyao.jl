
# Benchmark functions
# -------------------
function cpu_benchmark_fun(n)
    reg = rand(2n, 2n) |> MajoranaReg
    covmat = covariance_matrix(reg)
end

function gpu_benchmark_fun(n)
    reg = CUDA.rand(2n, 2n) |> MajoranaReg
    covmat = CUDA.@sync covariance_matrix(reg)
end
#=
function cpu_benchmark_fun(n)
    reg = rand(2n, 2n) |> MajoranaReg
    i = rand(1:n)
    pi = rand()
    ni = rand(Bool)
    #covmat = covariance_matrix(reg)
    update_covariance_matrix!(reg.state, i, pi, ni)
end

function gpu_benchmark_fun(n)
    reg = CUDA.rand(2n, 2n) |> MajoranaReg
    i = rand(1:n)
    pi = rand()
    ni = rand(Bool)
    #covmat = CUDA.@sync covariance_matrix(reg)
    CUDA.@sync FLOYao.update_covariance_matrix!(reg.state, i, pi, ni)
end
=#
#=
function update_covariance_matrix!(M::CuMatrix, i, pi, ni)
    n = size(M,1)
    @inline function kernel!(M)
        p, q = blockIdx().x + 2i, threadIdx().x + 2i + 1
        offset = (-1)^ni * (M[2i-1,q] * M[2i,p] - M[2i-1,p] * M[2i,q]) / (2pi)
        M[p,q] += ifelse(q>p, offset, zero(eltype(M)))
        return nothing
    end
    @cuda threads=n-2i-1 blocks=n-2i kernel!(M)
    return M
end
=#


