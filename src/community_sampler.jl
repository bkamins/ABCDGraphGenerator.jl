using Random
using StatsBase
using NearestNeighbors

function sample_points(n, d)
    points = randn(n, d)
    points ./= sqrt.(sum(x -> x^2, points, dims=2))
    points .*= rand(n) .^ (1/d)
    return (points .+ 1.0) ./ 2.0 # put the circle in the [0, 1]^d cube
end

# performance optimization future code
# function partition_points(x)
#     n, d = size(x)
#     cell_count = n / 10
#     dim_size = cell_count ^ (1/d)
#     side_mul = 2.0 ^ floor(log(dim_size) / log(2))
#     x_cut = floor.(Int32, x .* side_mul) .+ Int32(1)
#     x_cut_tuple = Tuple.(eachrow(x_cut))
#     x_partition = Dict{eltype(x_cut_tuple), Vector{Int32}}()
#     for (i, v) in enumerate(x_cut_tuple)
#         v = get!(x_partition, v) do
#             Int32[]
#         end
#         push!(v, i)
#     end
#     return x_partition, side_mul
# end

# get_coordinate(x) = floor.(Int32, x .* side_mul) .+ Int32(1)

function assign_points(x, c, p)
    @assert ndims(x) == 2
    @assert sum(c) == size(x, 1)
    @assert length(c) == length(p)
    x = copy(x)
    @assert size(x, 1) < typemax(Int32)
    all_idxs = collect(Int32, 1:size(x, 1))
    dist = vec(sum(x -> x^2, x, dims=2))
    res = Vector{Vector{Int32}}(UndefInitializer(), length(c))
    to_keep = trues(size(x, 1))
    for idx in p
        com = c[idx]
        ind = argmax(dist)
        ref = x[ind:ind, :]
        dist_c = vec(sum((x .- ref) .^ 2; dims=2))
        idxs = partialsortperm(dist_c, 1:com)
        res[idx] = all_idxs[idxs]
        to_keep[idxs] .= false
        x = x[to_keep, :]
        dist = dist[to_keep]
        all_idxs = all_idxs[to_keep]
        deleteat!(to_keep, sort!(idxs))
    end
    @assert size(x, 1) == 0
    @assert length(all_idxs) == 0
    @assert length(dist) == 0
    # TODO: disable below for production as overly expensive
    @assert sort(union(res...)) == 1:sum(c)
    return res
end

function assign_points2(kdtree, x, c, p)
    @assert ndims(x) == 2
    @assert sum(c) == size(x, 1)
    @assert length(c) == length(p)
    @assert size(x, 1) < typemax(Int32)

    res = Vector{Vector{Int32}}(UndefInitializer(), length(c))
    visited = falses(size(x, 1))
    dist = vec(sum(x -> x^2, x, dims=2))

    was_visited(idx) = visited[idx]

    for idx in p
        com = c[idx]
        ind = argmax(dist)
        ref = x[ind, :]
        idxs, _ = knn(kdtree, ref, Int(com), false, was_visited)
        res[idx] = Int32.(idxs)
        visited[idxs] .= true
        dist[idxs] .= -1.0
    end
    # TODO: disable below for production as overly expensive
    #   @assert sort(union(res...)) == 1:sum(c)
    return res
end

# note that this function returns node numbers from 1 to number_of_non_outlier_nodes
# for each community we get a set of nodes assigned to it
function populate_overlapping_clusters(coms::Vector{Int32}, η::Float64, d::Int)
    @assert 1 <= d
    true_coms = coms[2:end] # we are interested only in non-outlier communities
    grow_coms = [randround(s * η) for s in true_coms] # this is a target size of communities, as coms is primary community sizes
    p = randperm(length(true_coms)) # order in which communities are handled
    n = sum(true_coms)
    x = sample_points(n, d)

    #x_part = partition_points(x)
#   a2 = assign_points(x, true_coms, p)
    x2 = transpose(x)
    kdtree = KDTree(x2)
    a = assign_points2(kdtree, x, true_coms, p)
#    @assert length(a) == length(a2)
    # for (p1, q1) in zip(a, a2)
    #     @assert sort(p1) == sort(q1)
    # end
    @assert length.(a) == true_coms
    @assert sum(length, a) == sum(true_coms)

    # below we grow communities
    @assert length(a) == length(grow_coms)
    for (com, target) in zip(a, grow_coms)
        community_center = vec(mean(x[com, :], dims=1))
        ordering, _ = knn(kdtree, community_center, min(2 * target, size(x, 1)), true)
        com_set = Set(com)
        loc = 1
        while length(com) < target
            if loc > length(ordering)
                if length(ordering) < length(distances)
                    ordering2 = sortperm(distances)
                    @assert ordering == ordering2[1:length(ordering)]
                    ordering = ordering2
                else
                    throw(ArgumentError("η was too large"))
                end
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