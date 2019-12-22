using ABCDGraphGenerator

@info "Usage: julia graph_sampler.jl networkfile communityfile degreefile communitysizesfile μ isCL islocal"
@info "Example: julia graph_sampler.jl network.dat community.dat degrees.dat community_sizes.dat 0.2 true true"

networkfile = ARGS[1]
communityfile = ARGS[2]
degreefile = ARGS[3]
communitysizesfile = ARGS[4]
μ = parse(Float64, ARGS[5])
isCL = parse(Bool, ARGS[6])
islocal = parse(Bool, ARGS[7])

p = ABCDGraphGenerator.ABCDParams(parse.(Int, readlines(degreefile)),
                                  parse.(Int, readlines(communitysizesfile)),
                                  μ, isCL, islocal)

edges, clusters = gen_benchmark(p)

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
