using Pkg
using ABCDGraphGenerator

@info "Usage: julia abcd_sampler.jl config_filename"
@info "For the syntax of config_filename see example_config.toml file"

filename = ARGS[1]
conf = Pkg.TOML.parsefile(filename)

n = parse(Int, conf["n"])

τ₁ = parse(Float64, conf["t1"])
d_min = parse(Int, conf["d_min"])
d_max = parse(Int, conf["d_max"])
d_max_iter = parse(Int, conf["d_max_iter"])
@info "Expected value of degree: $(ABCDGraphGenerator.get_ev(τ₁, d_min, d_max))"
degs = ABCDGraphGenerator.sample_degrees(τ₁, d_min, d_max, n, d_max_iter)
open(io -> foreach(d -> println(io, d), degs), conf["degreefile"], "w")

τ₂ = parse(Float64, conf["t2"])
c_min = parse(Int, conf["c_min"])
c_max = parse(Int, conf["c_max"])
c_max_iter = parse(Int, conf["c_max_iter"])
@info "Expected value of community size: $(ABCDGraphGenerator.get_ev(τ₂, c_min, c_max))"
coms = ABCDGraphGenerator.sample_communities(τ₂, c_min, c_max, n, max_iter)
open(io -> foreach(d -> println(io, d), coms), conf["communitysizesfile"], "w")

ξ = parse(Float64, "xi")
isCL = parse(Bool, conf["isCL"])
p = ABCDGraphGenerator.ABCDParams(degs, coms, nothing, ξ, isCL, false)
edges, clusters = ABCDGraphGenerator.gen_graph(p)
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
