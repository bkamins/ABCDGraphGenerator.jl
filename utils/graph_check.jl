using Statistics

@info "Example sage: julia graph_check.jl degrees.dat community_sizes.dat community.dat network.dat [isCL]"

degrees_fname = ARGS[1]
community_sizes_fname = ARGS[2]
community_fname = ARGS[3]
network_fname = ARGS[4]
isCL = parse(Bool, ARGS[5])

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

@info "mean graph level internal fraction: $(mean(internal_frac))"

@info "Internal fractions per community:"
for i in 1:length(community_sizes)
    internal_frac_com = sum(internal_count[community .== i]) ./ sum((internal_count .+ outside_count)[community .== i])
    @info "Community $i has size $(community_sizes[i]) and internal fraction $internal_frac_com"
end
