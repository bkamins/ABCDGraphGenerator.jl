using ABCDGraphGenerator
using Random

# note that for backward compatibility reasons `[nout]` is an optional parameter
# that comes last
@info "Usage: julia graph_sampler.jl networkfile communityfile degreefile communitysizesfile mu|xi fraction isCL islocal [seed] [nout]"
@info "Example: julia graph_sampler.jl network.dat community.dat degrees.dat community_sizes.dat xi 0.2 true true 42 100"

networkfile = ARGS[1]
communityfile = ARGS[2]
degreefile = ARGS[3]
communitysizesfile = ARGS[4]
muxi = ARGS[5]
fraction = parse(Float64, ARGS[6])
isCL = parse(Bool, ARGS[7])
islocal = parse(Bool, ARGS[8])

length(ARGS) >= 9 &&  Random.seed!(parse(Int, ARGS[9]))
if length(ARGS) >= 10
    nout = parse(Float64, ARGS[10])
end

length(ARGS) >= 11 && @warn "more than 10 parameters passed"

coms = parse.(Int, readlines(communitysizesfile))

if nout > 0
    nout == coms[1] || throw(ArgumentError("nout does not match first community"))
end

muxi in ["mu","xi"] || throw(ArgumentError("only mu or xi names are allowed for"))
μ, ξ = nothing, nothing
if muxi == "mu"
    μ = fraction
else
    ξ = fraction
end

if isnothing(ξ) && nout > 0
    throw(ArgumentError("μ is not supported with outliers"))
end

if islocal && nout > 0
    throw(ArgumentError("local graph is not supported with outliers"))
end

if isCL && nout > 0
    throw(ArgumentError("Chung-Lu graph is not supported with outliers"))
end

p = ABCDGraphGenerator.ABCDParams(parse.(Int, readlines(degreefile)),
                                  coms,
                                  μ, ξ, isCL, islocal, nout > 0)

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
