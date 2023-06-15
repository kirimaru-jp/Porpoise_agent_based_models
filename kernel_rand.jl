
print("\033c")


using CUDA
CUDA.versioninfo()
CUDA.allowscalar(false)

using BenchmarkTools

#region "Ready functions"

CUDA.@time CUDA.rand_logn(Float32, 1, 5; mean=2, stddev=20)
@benchmark CUDA.rand_logn(Float32, 1, 5; mean=2, stddev=20)


using BinomialGPU
CUDA.@time rand_binomial(3, count = 1470, prob = 0.34)

# https://github.com/JuliaGPU/BinomialGPU.jl/blob/5054f7c285d0729067757c0405bc31acd78629e5/test/runtests.jl#L204
A = CUDA.zeros(Int, 1024, 1024)
println("")
println("Benchmarking constant parameter array: should run in less than 2ms on an RTX20xx card")
n = 128
p = 0.5
display(@benchmark CUDA.@sync rand_binomial!($A, count = $n, prob = $p))

B = CUDA.zeros(Int, 1000000)
display(@benchmark CUDA.@sync rand_binomial!($B, count = $n, prob = $p))
# BenchmarkTools.Trial: 10000 samples with 1 evaluation.
#  Range (min … max):  226.500 μs … 586.800 μs  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     301.300 μs               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   292.729 μs ±  31.949 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
#
#                                 ▄█▅▁
#   ▁▁▂▃▆▇▄▃▂▂▂▂▃▃▂▂▁▁▁▁▂▁▂▂▂▂▂▂▃▆████▇▅▃▂▃▂▂▂▂▂▂▂▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁ ▂
#   226 μs           Histogram: frequency by time          372 μs <
#
#  Memory estimate: 1.23 KiB, allocs estimate: 34.


#endregion


using Distributions, Random
a = CuArray{Int}(undef, 1024);

broadcast!(a) do
    dist = Binomial()
    rng = Random.default_rng()
    rand(rng, dist)
end




nothing


#region "CUDA"

A = CUDA.rand(10)

a = CuArray([1,2,3,4])
a .+= 1

for i in eachindex(a)
    a[i] += 1
end
a[1]


function vadd(c, a, b)
    i = threadIdx().x
    c[i] = a[i] + b[i]
    return
end

a = CuArray(1:10)
b = CuArray(2:2:20)
c = similar(a)
@cuda threads=length(a) vadd(c, a, b)
c


## Hardware limitations
threadIdx().x
a = CuArray(1:100000)
b = CuArray(2:2:200000)
c = similar(a)
@cuda threads=length(a) vadd(c, a, b)


CUDA.attribute(device(), CUDA.DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK)
function vadd(c, a, b)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= length(a)
        c[i] = a[i] + b[i]
    end
    return
end
@cuda threads=1024 blocks=cld(length(a),1024) vadd(c, a, b)
c


#endregion
