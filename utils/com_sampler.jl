using ABCDGraphGenerator

@info "Usage: julia com_sampler.jl filename τ₂ c_min c_max max_iter"
@info "Example: julia com_sampler.jl community_sizes.dat 3 5 50 1000"

τ₂ = parse(Float64, ARGS[1])
c_min = parse(Int, ARGS[2])
c_max = parse(Int, ARGS[3])
max_iter = parse(Int, ARGS[4])

@info "Expected value of community size: $(get_ev(τ₁, c_min, c_max))"

coms = sample_degrees(τ₂, c_min, c_max, n, max_iter)

open(io -> foreach(d -> println(io, d), coms), filename, "w")
