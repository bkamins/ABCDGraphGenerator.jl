using ABCDGraphGenerator
using Random

@info "Usage: julia com_sampler.jl filename τ₂ c_min c_max n max_iter [seed]"
@info "Example: julia com_sampler.jl community_sizes.dat 2 50 1000 10000 1000 42"

filename = ARGS[1]
τ₂ = parse(Float64, ARGS[2])
c_min = parse(Int, ARGS[3])
c_max = parse(Int, ARGS[4])
n = parse(Int, ARGS[5])
max_iter = parse(Int, ARGS[6])
length(ARGS) == 7 && Random.seed!(parse(Int, ARGS[7]))

@info "Expected value of community size: $(ABCDGraphGenerator.get_ev(τ₂, c_min, c_max))"

coms = ABCDGraphGenerator.sample_communities(τ₂, c_min, c_max, n, max_iter)

open(io -> foreach(d -> println(io, d), coms), filename, "w")
