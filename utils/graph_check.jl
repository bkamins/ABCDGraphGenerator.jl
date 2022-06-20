using Statistics
using ABCDGraphGenerator: ArgParse

function parse_commandline()
    s = ArgParse.ArgParseSettings()

    ArgParse.@add_arg_table! s begin
        "degrees"
            help = "degrees file"
            required = true
        "community_size"
            help = "community sizes file"
            required = true
        "community"
            help = "community file"
            required = true
        "network"
            help = "network file"
            required = true
        "isCL"
            help = "pass true if graph is CL and false if CM"
            arg_type = Bool
            required = true
    end
    s.usage = "graph_check.jl [-h] degrees community_size community network isCL"
    return ArgParse.parse_args(s)
end

parsed_args = parse_commandline()

degrees_fname = parsed_args["degrees"]
community_sizes_fname = parsed_args["community_size"]
community_fname = parsed_args["community"]
network_fname = parsed_args["network"]
isCL = parsed_args["isCL"]

degrees = parse.(Int, readlines(degrees_fname))
community_sizes = parse.(Int, readlines(community_sizes_fname))
community = (x -> parse.(Int, x[2])).(split.(readlines(community_fname)))
network = (x -> parse.(Int, x)).(split.(readlines(network_fname)))

@assert length(degrees) == length(community) == sum(community_sizes)
@info "Number of nodes: $(length(degrees))"
@info "Number of communities: $(length(community_sizes))"

nei_community = [Int[] for _ in 1:length(degrees)]
for (a, b) in network
    push!(nei_community[a], community[b])
    push!(nei_community[b], community[a])
end

@info "mean required degree: $(mean(degrees))"
@info "min required degree: $(minimum(degrees))"
@info "max required degree: $(maximum(degrees))"

@info "mean generated degree: $(mean(length.(nei_community)))"
@info "min generated degree: $(minimum(length.(nei_community)))"
@info "max generated degree: $(maximum(length.(nei_community)))"

if !isCL
    bad_degree = [i for i in 1:length(degrees) if degrees[i] != length(nei_community[i])]

    if isempty(bad_degree)
        @info "all generated degrees are equal to required degrees"
    else
        bad_nodes = [(node=i, expected=degrees[i], actual=length(nei_community[i])) for i in bad_degree]
        @warn "Nodes with not matching degrees are $bad_nodes"
    end
end

for i in 1:length(community_sizes)
    wanted_size = community_sizes[i]
    actual_size = count(==(i), community)
    if wanted_size != actual_size
        @warn "For community $i actual size $actual_size is not equal to wanted size $wanted_size"
    end
end

internal_count = [count(==(community[i]), nei_community[i]) for i in 1:length(degrees)]
outside_count = [count(!=(community[i]), nei_community[i]) for i in 1:length(degrees)]

internal_frac = internal_count ./ (internal_count .+ outside_count)

@info "mean graph level proportion of internal edges: $(mean(internal_frac))"

@info "Proportion of internal edges per community:"
for i in 1:length(community_sizes)
    internal_frac_com = sum(internal_count[community .== i]) ./ sum((internal_count .+ outside_count)[community .== i])
    @info "Community $i has size $(community_sizes[i]) and internal fraction $internal_frac_com"
end
