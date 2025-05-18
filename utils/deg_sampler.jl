using ABCDGraphGenerator
using Random

@info "Usage: julia deg_sampler.jl filename τ₁ d_min d_max n max_iter [seed]"
@info "Example: julia deg_sampler.jl degrees.dat 3 5 50 10000 1000 42"

filename = ARGS[1]
τ₁ = parse(Float64, ARGS[2])
d_min = parse(Int, ARGS[3])
d_max = parse(Int, ARGS[4])
n = parse(Int, ARGS[5])
max_iter = parse(Int, ARGS[6])
length(ARGS) == 7 && Random.seed!(parse(Int, ARGS[7]))

@info "Expected value of degree: $(ABCDGraphGenerator.get_ev(τ₁, d_min, d_max))"

degs = ABCDGraphGenerator.sample_degrees(τ₁, d_min, d_max, n, max_iter)

open(io -> foreach(d -> println(io, d), degs), filename, "w")
