using Pkg
using ABCDGraphGenerator
using Random

@info "Usage: julia abcd_oo_sampler.jl config_filename"
@info "For the syntax of config_filename see example_config_oo.toml file"

filename = ARGS[1]
conf = Pkg.TOML.parsefile(filename)
isempty(conf["seed"]) || Random.seed!(parse(Int, conf["seed"]))

nout = parse(Int, conf["nout"])
if nout < 0
    throw(ArgumentError("nout cannot be negative"))
end

ξ = parse(Float64, conf["xi"])
η = parse(Float64, conf["eta"])
if η < 1
    throw(ArgumentError("eta must be at least 1"))
end

d = parse(Int, conf["d"])
if d < 1
    throw(ArgumentError("d must be at least 1"))
end

ρ = parse(Float64, conf["rho"])

n = parse(Int, conf["n"])
if n < 0
    throw(ArgumentError("n cannot be negative"))
end

if nout > n
    throw(ArgumentError("number of outliers cannot be larger than graph size"))
end

# in what follows n is number of non-outlier nodes
n = n - nout

τ₁ = parse(Float64, conf["t1"])
d_min = parse(Int, conf["d_min"])
d_max = parse(Int, conf["d_max"])
d_max_iter = parse(Int, conf["d_max_iter"])
@info "Expected value of degree: $(ABCDGraphGenerator.get_ev(τ₁, d_min, d_max))"
degs = ABCDGraphGenerator.sample_degrees_oo(τ₁, d_min, d_max, n + nout, d_max_iter)
open(io -> foreach(d -> println(io, d), degs), conf["degreefile"], "w")

τ₂ = parse(Float64, conf["t2"])
c_min = parse(Int, conf["c_min"])
c_max = parse(Int, conf["c_max"])
c_max_iter = parse(Int, conf["c_max_iter"])
@info "Expected value of community size: $(ABCDGraphGenerator.get_ev(τ₂, c_min, c_max))"
coms = ABCDGraphGenerator.sample_communities_oo(τ₂, ceil(Int, c_min / η), floor(Int, c_max / η), n, c_max_iter)
@assert sum(coms) == n
pushfirst!(coms, nout)

p = ABCDGraphGenerator.ABCDParamsOO(degs, coms, ξ, η, d, ρ)
edges, clusters = ABCDGraphGenerator.gen_graph_oo(p)
open(conf["networkfile"], "w") do io
    for (a, b) in sort!(collect(edges))
        println(io, a, "\t", b)
    end
end
open(conf["communityfile"], "w") do io
    for (i, c) in enumerate(clusters)
        println(io, i, "\t", c)
    end
end

open(conf["communitysizesfile"], "w") do io
    comm_count = zeros(Int, length(coms))
    for c in clusters
        @assert length(c) > 0
        if 1 in c
            @assert length(c) == 1
        end
        for v in c
            comm_count[v] += 1
        end
    end
    println("eta is $η and empirically we have scaling of: ", extrema(comm_count[2:end] ./ coms[2:end]))
    foreach(d -> println(io, d), comm_count)
end
