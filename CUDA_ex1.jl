"""

Simplified version of BinomialGPU.jl

"""

print("\033c")


using CUDA
using Random

using CUDA: cuda_rng, i32
using Test
using BenchmarkTools

## extend the CUDA.jl functionality (rand, randn, rand_poisson, etc.) to include binomial distributions

const BinomialType = Union{Type{<:Integer}}
const BinomialArray = AnyCuArray{<:Integer}


function stirling_approx_tail(k)::Float32
    if k == 0
        return 0.0810614667953272f0
    elseif k == 1
        return 0.0413406959554092f0
    elseif k == 2
        return 0.0276779256849983f0
    elseif k == 3
        return 0.0207906721037650f0
    elseif k == 4
        return 0.0166446911898211f0
    elseif k == 5
        return 0.0138761288230707f0
    elseif k == 6
        return 0.0118967099458917f0
    elseif k == 7
        return 0.0104112652619720f0
    elseif k == 8
        return 0.00925546218271273f0
    elseif k == 9
        return 0.00833056343336287f0
    end
    kp1sq = (k + 1f0)^2;
    return (1.0f0 / 12 - (1.0f0 / 360 - 1.0f0 / 1260 / kp1sq) / kp1sq) / (k + 1)
end


function kernel_BTRS_scalar!(A, n, p, seed::UInt32, counter::UInt32)
    device_rng = Random.default_rng()

    # initialize the state
    @inbounds Random.seed!(device_rng, seed, counter)

    # grid-stride loop
    tid    = threadIdx().x
    window = blockDim().x * gridDim().x
    offset = (blockIdx().x - 1i32) * blockDim().x

    k = 0
    while offset < length(A)
        i = tid + offset

        r       = p/(1f0-p)
        s       = p*(1f0-p)

        stddev  = sqrt(n * s)
        b       = 1.15f0 + 2.53f0 * stddev
        a       = -0.0873f0 + 0.0248f0 * b + 0.01f0 * p
        c       = n * p + 0.5f0
        v_r     = 0.92f0 - 4.2f0 / b

        alpha   = (2.83f0 + 5.1f0 / b) * stddev;
        m       = floor((n + 1) * p)

        ks = 0f0

        while true
            usample = rand(Float32) - 0.5f0
            vsample = rand(Float32)

            us = 0.5f0 - abs(usample)
            ks = floor((2 * a / us + b) * usample + c)

            if us >= 0.07f0 && vsample <= v_r
                break
            end

            if ks < 0 || ks > n
                continue
            end

            v2 = CUDA.log(vsample * alpha / (a / (us * us) + b))
            ub = (m + 0.5f0) * CUDA.log((m + 1) / (r * (n - m + 1))) +
                    (n + 1) * CUDA.log((n - m + 1) / (n - ks + 1)) +
                    (ks + 0.5f0) * CUDA.log(r * (n - ks + 1) / (ks + 1)) +
                    stirling_approx_tail(m) + stirling_approx_tail(n - m) - stirling_approx_tail(ks) - stirling_approx_tail(n - ks)
            if v2 <= ub
                break
            end
        end

        if i <= length(A)
            @inbounds A[i] = ks
        end
        offset += window
    end
    return nothing
end


## constant (scalar) parameters
# function rand_binom!(A::BinomialArray, count::Integer, prob::AbstractFloat)
function rand_binom!(A, count::Integer, prob::AbstractFloat)
    n = count
    p = prob

	rng = cuda_rng()

    invert = false

	if p <= 0.5
        # BTRS algorithm for n*p > 10
            kernel = @cuda launch=false kernel_BTRS_scalar!(
                A, n, Float32(p), rng.seed, rng.counter
            )
    elseif p > 0.5
        invert = true
        p = 1 - p
        # BTRS algorithm for n*p > 10
            kernel = @cuda launch=false kernel_BTRS_scalar!(
                A, n, Float32(p), rng.seed, rng.counter
            )
    end

    config  = launch_configuration(kernel.fun)
    threads = max(32, min(config.threads, length(A)))
    blocks  = min(config.blocks, cld(length(A), threads))
    kernel(A, n, Float32(p), rng.seed, rng.counter; threads=threads, blocks=blocks)

    new_counter = Int64(rng.counter) + length(A)
    overflow, remainder = fldmod(new_counter, typemax(UInt32))
    rng.seed += overflow     # XXX: is this OK?
    rng.counter = remainder

    if invert
        return n .- A
    else
        return A
    end
end

# function rand_binomial!(rng, A::BinomialArray; count, prob)
#     return rand_binom!(rng, A, count, prob)
# end
# rand_binomial!(A::BinomialArray; kwargs...) = rand_binomial!(cuda_rng(), A; kwargs...)


n = 128
p = 0.5
B = CUDA.zeros(Int, 1000000)

@test n > 17
@test 0 < p < 1
@test n * p >= 10f0

# rand_binomial!(cuda_rng(), B, count = n, prob = p)
rand_binom!(B, n, p)

# display(@benchmark CUDA.@sync rand_binomial!($B, count = $n, prob = $p))
display(@benchmark CUDA.@sync rand_binom!($B, $n, $p))
# BenchmarkTools.Trial: 10000 samples with 1 evaluation.
#  Range (min … max):  233.200 μs …  11.261 ms  ┊ GC (min … max): 0.00% … 76.42%
#  Time  (median):     324.000 μs               ┊ GC (median):    0.00%    
#  Time  (mean ± σ):   326.702 μs ± 113.217 μs  ┊ GC (mean ± σ):  0.26% ±  
# 0.76%
#
#    ▁▂▁▂▁   ▁▂▁             ▁▅▇▇██▇▆▆▅▄▄▃▃▃▂▂▂▁▁▁                ▂
#   ▇█████▇█████▇▅▆▅▂▅▂▅▄▅▆▆▇███████████████████████▇█▇▇██▇▇▇▇▆▇▆ █        
#   233 μs        Histogram: log(frequency) by time        416 μs <
#
#  Memory estimate: 432 bytes, allocs estimate: 9.

# [https://cuda.juliagpu.org/stable/usage/memory/]
Array(B)
