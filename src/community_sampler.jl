using Random
using StatsBase

function sample_points(n)
    points = randn(n, 2)
    points ./= sqrt.(sum(x -> x^2, points, dims=2))
    points .*= rand(n) .^ 0.5
    return points
end

function assign_points(x, c, p)
    @assert ndims(x) == 2
    @assert sum(c) == size(x, 1)
    @assert length(c) == length(p)
    x = copy(x)
    all_idxs = collect(1:size(x, 1))
    dist = vec(sum(x -> x^2, x, dims=2))
    res = Vector{Vector{Int}}(UndefInitializer(), length(c))
    for idx in p
        com = c[idx]
        ind = argmax(dist)
        ref = x[ind, :]
        dist_c = [sum(x -> x^2, (r - ref)) for r in eachrow(x)]
        idxs = partialsortperm(dist_c, 1:com)
        res[idx] = all_idxs[idxs]
        to_keep = setdiff(1:size(x, 1), idxs)
        x = x[to_keep, :]
        dist = dist[to_keep]
        all_idxs = all_idxs[to_keep]
    end
    @assert size(x, 1) == 0
    @assert length(all_idxs) == 0
    @assert length(dist) == 0
    @assert sort(union(res...)) == 1:sum(c)
    return res
end

# note that this function returns node numbers from 1 to number_of_non_outlier_nodes
# for each community we get a set of nodes assigned to it
function populate_overlapping_clusters(coms::Vector{Int}, η::Float64)
    true_coms = coms[2:end] # we are interested only in non-outlier communities
    grow_coms = [randround(s * η) for s in true_coms] # this is a target size of communities, as coms is primary community sizes
    p = randperm(length(true_coms)) # order in which communities are handled
    n = sum(true_coms)
    x = sample_points(n)
    a = assign_points(x, true_coms, p)
    @assert length.(a) == true_coms
    @assert sum(length, a) == sum(true_coms)

    # below we grow communities
    @assert length(a) == length(grow_coms)
    for (com, target) in zip(a, grow_coms)
        community_center = vec(mean(x[com, :], dims=1))
        distances = [sum((v .- community_center) .^ 2) for v in eachrow(x)]
        ordering = sortperm(distances)
        com_set = Set(com)
        loc = 1
        while length(com) < target
            if loc > length(ordering)
                throw(ArgumentError("η was too large"))
            end
            point = ordering[loc]
            if !(point in com_set)
                push!(com, point)
            end
            loc += 1
        end
    end
    @assert length.(a) == grow_coms
    @assert sum(length, a) == sum(grow_coms)
    @assert all(allunique(c) for c in a)
    return [Set(c) for c in a]
end