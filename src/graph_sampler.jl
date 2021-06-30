"""
    ABCDParams

A structure holding parameters for ABCD graph generator. Fields:
* w::Vector{Int32}:             a sorted in descending order list of vertex degrees
* s::Vector{Int32}:             a sorted in descending order list of cluster sizes
* μ::Union{Float64, Nothing}: mixing parameter
* ξ::Union{Float64, Nothing}: background graph fraction
* isCL::Bool:                 if `true` a Chung-Lu model is used, otherwise configuration model
* islocal::Bool:              if `true` mixing parameter restriction is cluster local, otherwise
                              it is only global

Exactly one of ξ and μ must be passed as `Float64`. Also if `ξ` is passed then
`islocal` must be `false`.

The base ABCD graph is generated when ξ is passed and `isCL` is set to `false`.
"""
struct ABCDParams
    w::Vector{Int32}
    s::Vector{Int32}
    μ::Union{Float64, Nothing}
    ξ::Union{Float64, Nothing}
    isCL::Bool
    islocal::Bool

    function ABCDParams(w, s, μ, ξ, isCL, islocal)
        length(w) == sum(s) || throw(ArgumentError("inconsistent data"))
        if !isnothing(μ)
            0 ≤ μ ≤ 1 || throw(ArgumentError("inconsistent data on μ"))
        end
        if !isnothing(ξ)
            0 ≤ ξ ≤ 1 || throw(ArgumentError("inconsistent data ξ"))
            if islocal
                throw(ArgumentError("when ξ is provided local model is not allowed"))
            end
        end
        if isnothing(μ) && isnothing(ξ)
            throw(ArgumentError("inconsistent data: either μ or ξ must be provided"))
        end

        if !(isnothing(μ) || isnothing(ξ))
            throw(ArgumentError("inconsistent data: only μ or ξ may be provided"))
        end

        new(sort(w, rev=true),
            sort(s, rev=true),
            μ, ξ, isCL, islocal)
    end
end

function randround(x)
    d = floor(Int32, x)
    d + (rand() < x - d)
end

function populate_clusters(params::ABCDParams)
    w, s = params.w, params.s
    if isnothing(params.ξ)
        mul = 1.0 - params.μ
    else
        n = length(w)
        ϕ = 1.0 - sum((sl/n)^2 for sl in s)
        mul = 1.0 - params.ξ*ϕ
    end
    @assert length(w) == sum(s)
    @assert 0 ≤ mul ≤ 1
    @assert issorted(w, rev=true)
    @assert issorted(s, rev=true)

    slots = copy(s)
    clusters = Int32[]
    j = 0
    for (i, vw) in enumerate(w)
        while j + 1 ≤ length(s) && mul * vw + 1 ≤ s[j + 1]
            j += 1
        end
        j == 0 && throw(ArgumentError("could not find a large enough cluster for vertex of weight $vw"))
        wts = Weights(view(slots, 1:j))
        wts.sum == 0 && throw(ArgumentError("could not find an empty slot for vertex of weight $vw"))
        loc = sample(1:j, wts)
        push!(clusters, loc)
        slots[loc] -= 1
    end
    clusters
end

function CL_model(clusters, params)
    @assert params.isCL
    w, s, μ = params.w, params.s, params.μ
    cluster_weight = zeros(Int32, length(s))
    for i in axes(w, 1)
        cluster_weight[clusters[i]] += w[i]
    end
    total_weight = sum(cluster_weight)
    if params.islocal
        ξl = @. μ / (1.0 - cluster_weight / total_weight)
        maximum(ξl) >= 1 && throw(ArgumentError("μ is too large to generate a graph"))
    else
        if isnothing(params.ξ)
            ξg = μ / (1.0 - sum(x -> x^2, cluster_weight) / total_weight^2)
            ξg >= 1 && throw(ArgumentError("μ is too large to generate a graph"))
        else
            ξg = params.ξ
        end
    end

    wf = float.(w)
    edges = Set{Tuple{Int32, Int32}}()
    for i in axes(s, 1)
        local_edges = Set{Tuple{Int32, Int32}}()
        idxᵢ = findall(==(i), clusters)
        wᵢ = wf[idxᵢ]
        ξ = params.islocal ? ξl[i] : ξg
        m = randround((1-ξ) * sum(wᵢ) / 2)
        ww = Weights(wᵢ)
        while length(local_edges) < m
            a = sample(idxᵢ, ww, m - length(local_edges))
            b = sample(idxᵢ, ww, m - length(local_edges))
            for (p, q) in zip(a, b)
                p != q && push!(local_edges, minmax(p, q))
            end
        end
        union!(edges, local_edges)
    end
    wwt = if params.islocal
        Weights([ξl[clusters[i]]*x for (i,x) in enumerate(wf)])
    else
        Weights(ξg * wf)
    end
    while 2*length(edges) < total_weight
        a = sample(axes(w, 1), wwt, randround(total_weight / 2) - length(edges))
        b = sample(axes(w, 1), wwt, randround(total_weight / 2) - length(edges))
        for (p, q) in zip(a, b)
            p != q && push!(edges, minmax(p, q))
        end
    end
    edges
end

function config_model(clusters, params)
    @assert !params.isCL
    w, s, μ = params.w, params.s, params.μ

    cluster_weight = zeros(Int32, length(s))
    for i in axes(w, 1)
        cluster_weight[clusters[i]] += w[i]
    end
    total_weight = sum(cluster_weight)
    if params.islocal
        ξl = @. μ / (1.0 - cluster_weight / total_weight)
        maximum(ξl) >= 1 && throw(ArgumentError("μ is too large to generate a graph"))
        w_internal_raw = [w[i] * (1 - ξl[clusters[i]]) for i in axes(w, 1)]
    else
        if isnothing(params. ξ)
            ξg = μ / (1.0 - sum(x -> x^2, cluster_weight) / total_weight^2)
            ξg >= 1 && throw(ArgumentError("μ is too large to generate a graph"))
        else
            ξg = params.ξ
        end
        w_internal_raw = [w[i] * (1 - ξg) for i in axes(w, 1)]
    end

    clusterlist = [Int32[] for i in axes(s, 1)]
    for i in axes(clusters, 1)
        push!(clusterlist[clusters[i]], i)
    end
    # order by cluster size
    idx = sortperm([length(cluster) for cluster in clusterlist], rev=false)
    clusterlist = clusterlist[idx]

    edges::Vector{Set{Tuple{Int32, Int32}}} = []

    unresolved_collisions = 0
    w_internal = zeros(Int32, length(w_internal_raw))
    mutex = ReentrantLock()
    @threads for tid in 1:nthreads()
      local thr_clusters::Vector{Vector{Int32}} = []
      local thr_weights::Vector{Vector{Int32}} = []

      for c in tid:nthreads():length(s)
        local cluster = clusterlist[c]
        local w_cluster = zeros(Int32, length(cluster))
        local maxw_idx = argmax(view(w_internal_raw, cluster))
        local wsum = 0
        for i in axes(cluster, 1)
            if i != maxw_idx
                neww = randround(w_internal_raw[cluster[i]])
                w_cluster[i] = neww
                wsum += neww
            end
        end
        local maxw = floor(Int32, w_internal_raw[cluster[maxw_idx]])
        w_cluster[maxw_idx] = maxw + (isodd(wsum) ? iseven(maxw) : isodd(maxw))

        push!(thr_clusters, cluster)
        push!(thr_weights, w_cluster)
      end
      @debug "tid $(tid) getting 1 lock"
      lock(mutex)
        @debug "tid $(tid) got 1 lock"
        foreach((cluster,w_cluster)->w_internal[cluster]=w_cluster, thr_clusters, thr_weights)
        @debug "tid $(tid) releasing 1 lock"
      unlock(mutex)
    end

    @debug "GLOBAL starting"
    global_edges = Set{Tuple{Int32, Int32}}()
    recycle = Tuple{Int32,Int32}[]
    w_global = w - w_internal
    stubs::Vector{Int32} = zeros(Int32, sum(w_global))
    gt = Threads.@spawn begin
        sizehint!(global_edges, length(stubs)>>1)

        v::Vector{Int32} = cumsum(w_global)
        foreach((i,j,k)->stubs[i:j].=k, [1;v.+1], v, axes(w,1))
        @assert sum(w) == length(stubs) + sum(w_internal)
        shuffle!(stubs)

        @debug "$(length(edges)) communities"
        for i in 1:2:length(stubs)
            e = minmax(stubs[i], stubs[i+1])
            if (e[1] == e[2]) || (e in global_edges)
                push!(recycle, e)
            else
                push!(global_edges, e)
            end
        end
        @debug "dups1 are $(recycle)"
    end
    length_recycle = length(recycle)

    @threads for tid in 1:max(1, nthreads()-1)
      local thr_edges   = Set{Tuple{Int32, Int32}}[]
      local thr_recycle = Vector{Tuple{Int32,Int32}}[]

      for c in tid:max(1, nthreads()-1):length(s)
        local cluster = clusterlist[c]
        local w_cluster = w_internal[cluster]

        @debug "tid $(tid) cluster $(length(cluster)) w_cluster $(sum(w_cluster))"
        local v::Vector{Int32} = cumsum(w_cluster)
        local stubs::Vector{Int32} = zeros(Int32, sum(w_cluster))
        foreach((i,j,k)->stubs[i:j].=k, [1;v.+1], v, cluster)
        @assert sum(w_cluster) == length(stubs)

        shuffle!(stubs)
        local local_edges = Set{Tuple{Int32, Int32}}()
        sizehint!(local_edges, length(stubs)>>1)
        local recycle = Tuple{Int32,Int32}[]
        for i in 1:2:length(stubs)
            e = minmax(stubs[i], stubs[i+1])
            if (e[1] == e[2]) || (e in local_edges)
                push!(recycle, e)
            else
                push!(local_edges, e)
            end
        end
        local last_recycle = length(recycle)
        local recycle_counter = last_recycle
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
            local p1 = popfirst!(recycle)
            local from_recycle = 2 * length(recycle) / length(stubs)
            local success = false
            for _ in 1:2:length(stubs)
                local p2 = if rand() < from_recycle
                    used_recycle = true
                    recycle_idx = rand(axes(recycle, 1))
                    recycle[recycle_idx]
                else
                    used_recycle = false
                    rand(local_edges)
                end
                if rand() < 0.5
                    local newp1 = minmax(p1[1], p2[1])
                    local newp2 = minmax(p1[2], p2[2])
                else
                    local newp1 = minmax(p1[1], p2[2])
                    local newp2 = minmax(p1[2], p2[1])
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
        push!(thr_edges, local_edges)
        push!(thr_recycle, recycle)
      end
      @debug "tid $(tid) getting 2 lock"
      lock(mutex)
        @debug "tid $(tid) got 2 lock"
        append!(edges, thr_edges)
        append!(recycle, thr_recycle...)
        @debug "tid $(tid) releasing 2 lock"
      unlock(mutex)
    end

    unresolved_collisions = length(recycle) - length_recycle
    if unresolved_collisions > 0
        println("Unresolved_collisions: ", unresolved_collisions,
                "; fraction: ", 2 * unresolved_collisions / total_weight)
    end

    @debug "GLOBAL waiting to complete"
    wait(gt)
    @debug "GLOBAL resolving dups"
    @debug "intersect $(length(global_edges)) global_edges with $(typeof(edges)) $([length(e) for e in edges])"
    dups = [Set{Tuple{Int32, Int32}}() for _ in axes(edges, 1)]
    @threads for i in axes(edges, 1)
        if length(global_edges) > length(edges[i])
            dups[i] = intersect(global_edges, edges[i])
        else
            dups[i] = intersect(edges[i], global_edges)
        end
    end
    append!(recycle, dups...)
    setdiff!(global_edges, dups...)
    dups = Nothing
    @debug "dups2 are $(recycle)"
    while !isempty(recycle)
        p1 = pop!(recycle)
        from_recycle = 2 * length(recycle) / length(stubs)
        p2 = if rand() < from_recycle
            i = rand(axes(recycle, 1))
            recycle[i], recycle[end] = recycle[end], recycle[i]
            pop!(recycle)
        else
            x = rand(global_edges)
            pop!(global_edges, x)
        end
        if rand() < 0.5
            newp1 = minmax(p1[1], p2[1])
            newp2 = minmax(p1[2], p2[2])
        else
            newp1 = minmax(p1[1], p2[2])
            newp2 = minmax(p1[2], p2[1])
        end
        for newp in (newp1, newp2)
            if (newp[1] == newp[2]) || (newp in global_edges) || any(cluster->newp in cluster, edges)
                push!(recycle, newp)
            else
                push!(global_edges, newp)
            end
        end
    end
    @debug "dups3 are $(recycle)" # should be empty
    # old_len = length(edges)
    push!(edges, global_edges)
    # @assert length(edges) == old_len + length(global_edges)
    @debug "$(length(global_edges)) global_edges $(length(stubs)) stubs"
    @assert 2 * length(global_edges) == length(stubs)
    ChainedVector([edgeset.dict.keys[edgeset.dict.slots.==0x1] for edgeset in edges])
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
    edges = params.isCL ? CL_model(clusters, params) : config_model(clusters, params)
    (edges=edges, clusters=clusters)
end
