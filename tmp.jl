using CUDA, BenchmarkTools

function swap_and_sign!(A::Matrix)
    n = size(A, 1) ÷ 2
    for i in 1:2n
        for j in 2:2:2n
            A[i,j], A[i,j-1] = A[i,j-1], -A[i,j]
        end
    end
    nothing
end

function swap_and_sign!(A::CuArray)
    n = size(A, 1) ÷ 2
    threads = n
    blocks = 2n
    @inline function kernel(data)
        i, j = blockIdx().x, 2 * threadIdx().x
        data[i,j], data[i,j-1] = data[i,j-1], -data[i,j]
        nothing
    end
    @cuda threads=threads blocks=blocks kernel(A)
    nothing
end

function swap_and_sign_views!(A::AbstractMatrix)
    n = size(A, 1) ÷ 2
    for j in 2:2:2n
        tmp = A[:,j]
        A[:,j] .= A[:,j-1]
        A[:,j-1] .= -tmp
    end
    nothing
end

arr = rand(1000, 1000)
cuarr = arr |> cu


@btime CUDA.@sync swap_and_sign!(cuarr)
#  2.469 ms (19 allocations: 576 bytes)

@btime swap_and_sign!(arr)
#  1.563 ms (0 allocations: 0 bytes)
