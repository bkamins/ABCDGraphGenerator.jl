using ABCDGraphGenerator
using Random

@info "Usage: julia graph_sampler.jl networkfile communityfile degreefile communitysizesfile mu|xi fraction isCL islocal [seed]"
@info "Example: julia graph_sampler.jl network.dat community.dat degrees.dat community_sizes.dat xi 0.2 true true 42"

networkfile = ARGS[1]
communityfile = ARGS[2]
degreefile = ARGS[3]
communitysizesfile = ARGS[4]
muxi = ARGS[5]
fraction = parse(Float64, ARGS[6])
isCL = parse(Bool, ARGS[7])
islocal = parse(Bool, ARGS[8])

muxi in ["mu","xi"] || throw(ArgumentError("only mu or xi names are allowed for"))
μ, ξ = nothing, nothing
if muxi == "mu"
    μ = fraction
else
    ξ = fraction
end
length(ARGS) == 9 &&  Random.seed!(parse(Int, ARGS[9]))

p = ABCDGraphGenerator.ABCDParams(parse.(Int, readlines(degreefile)),
                                  parse.(Int, readlines(communitysizesfile)),
                                  μ, ξ, isCL, islocal)

edges, clusters = ABCDGraphGenerator.gen_graph(p)

open(networkfile, "w") do io
    for (a, b) in sort!(collect(edges))
        println(io, a, "\t", b)
    end
end

open(communityfile, "w") do io
    for (i, c) in enumerate(clusters)
        println(io, i, "\t", c)
    end
end
