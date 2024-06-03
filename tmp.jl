using CUDA, BenchmarkTools
using KernelAbstractions

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
    @inline function kernel_fun!(data)
        i = (blockIdx().x-1) * blockDim().x + threadIdx().x
        i_stride = gridDim().x * blockDim().x
        
        # multiply with 2 here for correct indexing.
        j = 2*((blockIdx().y-1) * blockDim().y + threadIdx().y) 
        j_stride = 2 * gridDim().y * blockDim().y

        while i <= 2n
            while j <= 2n
                data[i,j], data[i,j-1] = data[i,j-1], -data[i,j]
                j += j_stride
            end
            i += i_stride
        end
        nothing
    end
    kernel = @cuda launch=false kernel_fun!(A)
    config = launch_configuration(kernel.fun)

    # this seemed to give the faster iteration order
    threads_x = min(2n, config.threads)
    threads_y = min(fld(config.threads, threads_x), n) # put all remaining threads in here
    blocks_y = cld(n, threads_y)
    blocks_x = cld(2n, threads_x)
    threads = (threads_x, threads_y)
    blocks = (blocks_x, blocks_y)
    CUDA.@sync kernel(A; threads, blocks)
    nothing
end

function swap_and_sign_KA!(A::CuArray)
    backend = KernelAbstractions.get_backend(A)
    kernel! = swap_and_sign_kernel!(backend)
    kernel!(A, ndrange=(size(A, 1) , size(A,2) ÷ 2))
end

@kernel function swap_and_sign_kernel!(A)
    i, j = @index(Global, NTuple)
    j = 2 * j
    A[i,j], A[i,j-1] = A[i,j-1], -A[i,j]
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
