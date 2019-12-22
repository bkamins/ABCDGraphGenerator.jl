using ABCDGraphGenerator

@info "Usage: julia deg_sampler.jl filename τ₁ d_min d_max max_iter"
@info "Example: julia deg_sampler.jl degrees.dat 3 5 50 1000"

τ₁ = parse(Float64, ARGS[1])
d_min = parse(Int, ARGS[2])
d_max = parse(Int, ARGS[3])
max_iter = parse(Int, ARGS[4])

@info "Expected value of degree: $(get_ev(τ₁, d_min, d_max))"

degs = sample_degrees(τ₁, d_min, d_max, n, max_iter)

open(io -> foreach(d -> println(io, d), degs), filename, "w")
