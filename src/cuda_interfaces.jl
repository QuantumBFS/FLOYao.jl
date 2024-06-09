"""
    FLOYao.cu(reg)

Upload the register state from CPU to GPU.
"""
function cu end

"""
    FLOYao.cpu(cureg)

Download the register state from GPU to CPU.
"""
function cpu end

"""
    FLOYao.cuzero_state([T=Float32,] n)

The GPU version of [`FLOYao.zero_state`](@ref).
"""
function cuzero_state end

"""
    FLOYao.cuone_state([T=Float32,] n)

The GPU version of [`FLOYao.one_state`](@ref).
"""
function cuone_state end

"""
    FLOYao.cuproduct_state([T=Float64,] bit_str::DitStr{2})
    FLOYao.cuproduct_state([T=Float64,] bit_configs::AbstractVector)
    FLOYao.cuproduct_state([T=Float64,] nbits::Int, val::Int)

The GPU version of [`FLOYao.product_state`](@ref)
"""
function cuproduct_state end
