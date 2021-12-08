using StatsBase, Random

function trunc_powerlaw_weigths(α::Real, v_min::Int, v_max::Int)
    @assert α >= 1
    @assert 1 <= v_min <= v_max
    w = Weights([1/i^α for i in v_min:v_max])
end

function sample_trunc_powerlaw(α::Real, v_min::Int, v_max::Int, n::Int)
    @assert α >= 1
    @assert 1 <= v_min <= v_max
    @assert n > 0
    w = Weights([1/i^α for i in v_min:v_max])
    sample(v_min:v_max, w, n)
end

sample_trunc_powerlaw(W::Weights, v_min::Int, v_max::Int, n::Int) = (@assert n > 0; sample(v_min:v_max, W, n))

"""
    get_ev(α, v_min, v_max)

Return the expected value of truncated discrete power law distribution
with truncation range `[v_min, v_max]` and exponent `α`.
"""
function get_ev(α, v_min, v_max)
    @assert α > 1
    @assert 1 <= v_min <= v_max
    w = [1/i^α for i in v_min:v_max]
    w /= sum(w)
    sum(x*w[i] for (i, x) in enumerate(v_min:v_max))
end

"""
    find_v_min(α, v_max, v_ev)

Return the lower truncation given expected value `v_ev` and upper truncation `v_max`
of truncated discrete power law distribution with exponent `α`.
"""
function find_v_min(α, v_max, v_ev)
    @assert v_ev <= v_max
    v_min = v_ev
    cur_err = abs(v_ev, get_ev(α, v_min, v_max))
    while true
        v_min == 1 && return 1
        v_min -= 1
        next_err = abs(v_ev, get_ev(α, v_min, v_max))
        if next_err > cur_err
            return v_min + 1
        else
            cur_err = next_err
        end
    end
end

"""
    sample_degrees(τ₁, d_min, d_max, n, max_iter)

Return a vector of length `n` of sampled degrees of vertices following a truncated
discrete power law distribution with truncation range `[d_min, d_max]` and exponent `τ₁`.

The sampling is attempted `max_iter` times to find an admissible result.
In case of failure a correction of degree is applied to a degree of the vertex
with highest degree.

The producedure does not check if the returned vector is a graphical degree sequence.
"""
function sample_degrees(τ₁, d_min, d_max, n, max_iter)
    local s
    for i in 1:max_iter
        s = sample_trunc_powerlaw(τ₁, d_min, d_max, n)
        iseven(sum(s)) && return sort!(s, rev=true)
    end
    @warn "Failed to sample an admissible degree sequence in $max_iter draws. Fixing"
    i = argmax(s)
    if s[i] == 0
        s[i] = 1 # this should not happen in practice
    else
        s[i] -= 1
    end
    return sort!(s, rev=true)
end

"""
    sample_communities(τ₂, c_min, c_max, n, max_iter)

Return a vector of sampled community sizes following a truncated
discrete power law distribution with truncation range `[c_min, c_max]` and exponent `τ₂`.
The sum of sizes is equal to number of vertices in the graph `n`.

The sampling is attempted `max_iter` times to find an admissible result.
In case of failure a correction of community sizes is applied to the sampled sequence
that was closest to a feasible one to ensure that the result is admissible.
"""
function sample_communities(τ₂, c_min, c_max, n, max_iter)
    @assert 1 <= c_min <= c_max
    l_min = n / c_max
    l_max = n / c_min
    @assert l_min >= 1
    @assert ceil(l_min) <= floor(l_max)
    local best_s
    local best_ss = typemax(Int)
    w = trunc_powerlaw_weigths(τ₂, c_min, c_max)
    for i in 1:max_iter
        s = sample_trunc_powerlaw(w, c_min, c_max, ceil(Int, l_max))
        stopidx = 0
        ss = 0
        while ss < n
            stopidx += 1
            ss += s[stopidx]
        end
        ss == n && return sort!(s[1:stopidx], rev=true)
        if ss < best_ss
            best_ss = ss
            best_s = s[1:stopidx]
        end
    end
    @warn "failed to sample an admissible community sequence in $max_iter draws. Fixing"
    shuffle!(best_s)
    if length(best_s) > l_max
        best_s = best_s[1:l_max]
        best_ss = sum(best_s)
    end
    i = 0
    while best_ss != n
        if i ≥ length(best_s)
            i = 0
            shuffle!(best_s)
        end
        i += 1
        change = sign(n - best_ss)
        if change > 0
            best_s[i] < c_max || continue
        else
            best_s[i] > c_min || continue
        end
        best_ss += change
        best_s[i] += change
    end
    return sort!(best_s, rev=true)
end
