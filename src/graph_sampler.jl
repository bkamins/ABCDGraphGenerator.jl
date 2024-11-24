"""
    ABCDParams

A structure holding parameters for ABCD graph generator. Fields:
* w::Vector{Int}:   a sorted in descending order list of vertex degrees
* s::Vector{Int}:   a sorted in descending order list of cluster sizes
                    (except first community which is outliers);
                    cluster sizes must be for primary community
* ξ::Float64:       background graph fraction
* η::Float64:       average number of communities a non-outlier node is part of; default 1

The graph will be generated using configuration model global approach.
"""
struct ABCDParams
    w::Vector{Int}
    s::Vector{Int}
    ξ::Float64
    η::Float64

    function ABCDParams(w::Vector{Int}, s::Vector{Int}, ξ::Float64, η::Float64)
        all(>(0), w) || throw(ArgumentError("all degrees must be positive"))
        length(w) == sum(s) || throw(ArgumentError("inconsistent data"))
        length(s) < 2 && throw(ArgumentError("no communities requested"))
        s[1] >= 0 || throw(ArgumentError("negative count of outliers passed"))
        0 ≤ ξ ≤ 1 || throw(ArgumentError("inconsistent data ξ"))
        η < 1 && throw(ArgumentError("η must be greater or equal than 1"))

        news = copy(s)
        all(>(0), @view(news[2:end])) || throw(ArgumentError("all community sizes must be positive"))
        sort!(@view(news[2:end]), rev=true)

        largest = news[2] # size of largest non-outlier community
        if η * largest > length(w) - news[1]
            throw(ArgumentError("η must be small enough so that overlapping communities are not too big"))
        end

        new(sort(w, rev=true), news, ξ, η)
    end
end

function randround(x)
    d = floor(Int, x)
    d + (rand() < x - d)
end

function populate_clusters(params::ABCDParams)
    w, s = params.w, params.s

    n = length(w)
    s0 = s[1]
    ϕ = 1.0 - sum((sl/(n-s0))^2 for sl in s[2:end]) * (n-s0)*params.ξ / ((n-s0)*params.ξ + s0)
    mul = 1.0 - params.ξ*ϕ

    @assert length(w) == sum(s)
    @assert 0 ≤ mul ≤ 1
    @assert issorted(w, rev=true)
    @assert issorted(s[2:end], rev=true)

    slots = copy(s) # number of slots left in a community to be assigned
    clusters = [Int[] for i in 1:length(w)] # primary cluster of a node, [1] is outlier community, Int[] is no community yet

    # handle outliers
    nout = s[1]
    n = length(params.w)
    L = sum(d -> min(1.0, params.ξ * d), params.w)
    threshold = L + nout - L * nout / n - 1.0 # we cannot put too heavy nodes as outliers
    idx = findfirst(<=(threshold), params.w)
    @assert all(i -> params.w[i] <= threshold, idx:n)
    if length(idx:n) < nout
        throw(ArgumentError("not enough nodes feasible for classification as outliers"))
    end
    tabu = sample(idx:n, nout, replace=false)
    for i in tabu
        push!(clusters[i], 1) # outlier community
    end
    stabu = Set(tabu) # stabu is a set of indices already used up

    # handle normal communities
    # note that numbers assigned to communities are from 1 to sum(slots[2:end]) so remapping is needed later
    slots_less_1 = slots[2:end]
    cluster_assignments = populate_overlapping_clusters(slots, params.η)
    ηu = [count(c -> i in c, cluster_assignments) for i in 1:sum(slots_less_1)] # count in how many communities each number belongs
    @assert minimum(ηu) >= 1
    min_com = [minimum(length(c) - 1 for c in cluster_assignments if i in c) for i in 1:sum(slots_less_1)] # minimum size of a community less one for each number
    max_degree = (ηu .* min_com) / mul

    for (i, vw) in enumerate(w)
        i in stabu && continue # skip nodes in outlier community
        good_idxs = findall(md -> vw <= md, max_degree) # later make it faster, but for now leave a simple implementation
        isempty(good_idxs) && throw(ArgumentError("could not find a large enough cluster for vertex of weight $vw with index $i"))
        chosen_idx = rand(good_idxs)
        clusters[i] = findall(c -> chosen_idx in c, cluster_assignments) .+ 1 # write down cluster numbers of chosen node; need to add 1 as first cluster is for outliers
        max_degree[chosen_idx] = -1 # make sure we will not use chosen_idx later; note that this needs refactoring if the code is optimized for speed later
    end

    @assert sum(length, clusters) == s0 + sum(length, cluster_assignments)

    return clusters # which clusters a given node is assigned to
end

function generate_initial_graph(weights::Vector{Int})
    stubs = Int[]
    for i in 1:length(weights)
        for _ in 1:weights[i]
            push!(stubs, i)
        end
    end

    @assert sum(weights) == length(stubs)
    @assert iseven(length(stubs))

    shuffle!(stubs)

    local_edges = Set{Tuple{Int,Int}}()
    recycle = Tuple{Int,Int}[]

    for i in 1:2:length(stubs)
        e = minmax(stubs[i], stubs[i+1])
        if (e[1] == e[2]) || (e in local_edges)
            push!(recycle, e)
        else
            push!(local_edges, e)
        end
    end
    last_recycle = length(recycle)
    recycle_counter = last_recycle
    while !isempty(recycle)
        recycle_counter -= 1
        if recycle_counter < 0
            if length(recycle) < last_recycle
                last_recycle = length(recycle)
                recycle_counter = last_recycle
            else
                break
            end
        end
        p1 = popfirst!(recycle)
        from_recycle = 2 * length(recycle) / length(stubs)
        success = false
        for _ in 1:2:length(stubs)
            p2 = if rand() < from_recycle
                used_recycle = true
                recycle_idx = rand(axes(recycle, 1))
                recycle[recycle_idx]
            else
                used_recycle = false
                rand(local_edges)
            end
            if rand() < 0.5
                newp1 = minmax(p1[1], p2[1])
                newp2 = minmax(p1[2], p2[2])
            else
                newp1 = minmax(p1[1], p2[2])
                newp2 = minmax(p1[2], p2[1])
            end
            if newp1 == newp2
                good_choice = false
            elseif (newp1[1] == newp1[2]) || (newp1 in local_edges)
                good_choice = false
            elseif (newp2[1] == newp2[2]) || (newp2 in local_edges)
                good_choice = false
            else
                good_choice = true
            end
            if good_choice
                if used_recycle
                    recycle[recycle_idx], recycle[end] = recycle[end], recycle[recycle_idx]
                    pop!(recycle)
                else
                    pop!(local_edges, p2)
                end
                success = true
                push!(local_edges, newp1)
                push!(local_edges, newp2)
                break
            end
        end
        success || push!(recycle, p1)
    end

    unused_stubs = Int[]
    for (a, b) in recycle
        push!(unused_stubs, a, b)
    end

    @assert sum(weights) == length(local_edges) * 2 + length(unused_stubs)

    return local_edges, unused_stubs
end

function config_model(clusters, params)
    w, s, ξ = params.w, params.s, params.ξ

    @assert iseven(sum(w))
    w_internal_raw = randround.([w[i] * (1 - ξ) for i in axes(w, 1)])
    for i in findall(==([1]), clusters)
        w_internal_raw[i] = 0
    end

    w_external = w - w_internal_raw

    clusterlist = [Int[] for i in 1:maximum(c -> maximum(c), clusters)] # list of nodes in each cluster
    @show length(clusters)
    @show length(clusterlist)
    for i in axes(clusters, 1)
        c = clusters[i]
        for x in c
            push!(clusterlist[x], i)
        end
    end

    w_internal_comm = [zeros(Int, length(w_internal_raw)) for i in 1:length(clusterlist)] # this holds internal degree of each community

    for i in axes(clusters, 1)
        wi = w_internal_raw[i]
        nc = length(clusters[i])
        share = floor(Int, wi / nc)
        extra = wi - nc * share
        z = shuffle(1:nc)[1:extra]
        for j in 1:nc
            w_internal_comm[clusters[i][j]][i] = share + (j in z)
        end
    end

    for wic in w_internal_comm # make sure that for each community sum of its degrees is even
        if isodd(sum(wic))
            largest = argmax(wic)
            @assert wic[largest] > 0
            wic[largest] -= 1
            w_external[largest] += 1
        end
    end

    @assert sum(w_internal_comm) + w_external == w
    @assert iseven(length(w_external))
    @assert all(x -> iseven(length(x)), w_internal_comm)
    @assert all(==(0), w_internal_comm[1])


    partial_graphs = Set{Tuple{Int,Int}}[]
    unused_stubs = Int[]

    idxs_com = 0
    for w_int in w_internal_comm
        idxs_com += 1
        if idxs_com == 1 # outlier community
            @assert sum(w_int) == 0
        else
            g, s = generate_initial_graph(w_int)
            push!(partial_graphs, g)
            append!(unused_stubs, s)
        end
    end

    let
        g, s = generate_initial_graph(w_external)
        push!(partial_graphs, g)
        append!(unused_stubs, s)
    end

    edges = Set{Tuple{Int,Int}}()

    for g in partial_graphs
        for e in g
            if e in edges
                push!(unused_stubs, e[1], e[2]) # duplicate across subgraphs
            else
                push!(edges, e)
            end
        end
    end

    @assert sum(w) == length(edges) * 2 + length(unused_stubs)
    @assert iseven(length(unused_stubs))

    last_recycle_len = length(unused_stubs)
    if length(unused_stubs) > 0
        println("Duplicates or self loops generated in total number of $(div(length(unused_stubs), 2)). Fixing.")

        while true
            shuffle!(unused_stubs)
            recycle = Int[]
            for i in 1:2:length(unused_stubs)
                e = minmax(unused_stubs[i], unused_stubs[i+1])
                if (e[1] == e[2]) || (e in edges)
                    push!(recycle, e[1], e[2])
                else
                    push!(edges, e)
                end
            end
            if length(recycle) == last_recycle_len # no success in generating new feasible edges
                println("Could not fix creation of $(div(last_recycle_len, 2)) edges. Giving up.")
                break
            end
            last_recycle_len = length(recycle)
            unused_stubs = recycle
            @assert iseven(length(recycle))
            isempty(recycle) && break
        end
    end

    @assert sum(w) == length(edges) * 2 + last_recycle_len
    return edges
end

"""
    gen_graph(params::ABCDParams)

Generate ABCD graph following parameters specified in `params`.

Return a named tuple containing a set of edges of the graph and a list of cluster
assignments of the vertices.
The ordering of vertices and clusters is in descending order (as in `params`).
"""
function gen_graph(params::ABCDParams)
    clusters = populate_clusters(params)
    edges = config_model(clusters, params)
    (edges=edges, clusters=clusters)
end
