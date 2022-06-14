using ABCDGraphGenerator
using Random

# note that for backward compatibility reasons `[nout]` is an optional parameter
# that comes last
@info "Usage: julia com_sampler.jl filename τ₂ c_min c_max n max_iter [seed] [nout]"
@info "Example: julia com_sampler.jl community_sizes.dat 2 50 1000 10000 1000 42 40"

filename = ARGS[1]
τ₂ = parse(Float64, ARGS[2])
c_min = parse(Int, ARGS[3])
c_max = parse(Int, ARGS[4])
n = parse(Int, ARGS[5])
max_iter = parse(Int, ARGS[6])
length(ARGS) >= 7 && Random.seed!(parse(Int, ARGS[7]))
if length(ARGS) >= 8
    nout = parse(Int, ARGS[8])
else
    nout = 0
end

length(ARGS) >= 9 && @warn "more than 8 parameters passed"

@info "Expected value of community size: $(ABCDGraphGenerator.get_ev(τ₂, c_min, c_max))"

if nout > n
    throw(ArgumentError("number of outliers cannot be larger than graph size"))
end

coms = ABCDGraphGenerator.sample_communities(τ₂, c_min, c_max, n - nout, max_iter)

if nout > 0
    pushfirst!(coms, nout)
end

open(io -> foreach(d -> println(io, d), coms), filename, "w")
