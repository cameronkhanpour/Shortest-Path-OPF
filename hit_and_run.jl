module HitAndRun

using LinearAlgebra, Random, Printf

"""
    _random_direction(dim::Int, direction_basis)

Generates a random unit vector in `dim` dimensions. If `direction_basis` is
provided, the direction is sampled from the span of its columns; otherwise it is
sampled isotropically from a multivariate normal distribution.
"""
function _random_direction(dim::Int,
        direction_basis::Union{Nothing, Matrix{Float64}}=nothing)
    if isnothing(direction_basis)
        direction = randn(dim)
    else
        coeff = randn(size(direction_basis, 2))
        direction = direction_basis * coeff
        if norm(direction) < 1e-12
            direction = randn(dim)
        end
    end
    return direction / norm(direction)
end

"""
    _line_box_parameter_range(direction, point, box; tol)

Finds the parameter interval `[t_min, t_max]` for which the line
`point + t * direction` stays inside the hyper-rectangle `box`.
"""
function _line_box_parameter_range(direction::Vector{Float64},
        point::Vector{Float64},
        box::Vector{Tuple{Float64, Float64}};
        tol::Float64=1e-8)
    t_min = -Inf
    t_max = Inf

    for i in eachindex(direction)
        dir_i = direction[i]
        lo_i, hi_i = box[i]

        if abs(dir_i) <= tol
            if point[i] < lo_i - tol || point[i] > hi_i + tol
                return 0.0, 0.0, false
            end
            continue
        end

        t1 = (lo_i - point[i]) / dir_i
        t2 = (hi_i - point[i]) / dir_i
        t_lo = min(t1, t2)
        t_hi = max(t1, t2)
        t_min = max(t_min, t_lo)
        t_max = min(t_max, t_hi)

        if t_min > t_max + tol
            return 0.0, 0.0, false
        end
    end

    if !isfinite(t_min) || !isfinite(t_max) || (t_max - t_min) <= tol
        return 0.0, 0.0, false
    end

    return t_min, t_max, true
end

"""
    _locate_feasibility_transition(current_point, direction, t_left, left_flag, t_right, right_flag,
                                   is_feasible; tol, max_iter)

Refines a transition point between a feasible and an infeasible endpoint on the
same line segment.
"""
function _locate_feasibility_transition(current_point::Vector{Float64},
        direction::Vector{Float64},
        t_left::Float64,
        left_flag::Bool,
        t_right::Float64,
        right_flag::Bool,
        is_feasible::Function;
        tol::Float64=1e-8,
        max_iter::Int=60)
    left_flag == right_flag &&
        error("Transition refinement requires endpoints with different feasibility labels.")

    lo = t_left
    hi = t_right
    lo_flag = left_flag

    for _ in 1:max_iter
        mid = (lo + hi) / 2.0
        mid_flag = is_feasible(current_point .+ mid .* direction)
        if mid_flag == lo_flag
            lo = mid
        else
            hi = mid
        end
        abs(hi - lo) < tol && break
    end

    return (lo + hi) / 2.0
end

"""
    _feasible_intervals_along_line(current_point, direction, t_min, t_max, is_feasible; max_step, tol, max_iter)

Approximates the full set `S ∩ {current_point + t * direction : t ∈ [t_min, t_max]}`
as a union of feasible parameter intervals. The scan resolution is controlled by
`max_step`, so thin components or holes smaller than that scale may be missed.
"""
function _feasible_intervals_along_line(current_point::Vector{Float64},
        direction::Vector{Float64},
        t_min::Float64,
        t_max::Float64,
        is_feasible::Function;
        max_step::Float64=0.02,
        tol::Float64=1e-8,
        max_iter::Int=60)
    if (t_max - t_min) <= tol
        return Tuple{Float64, Float64}[]
    end

    n_steps = max(1, ceil(Int, (t_max - t_min) / max_step))
    grid = collect(range(t_min, t_max; length=n_steps + 1))
    if !any(abs(t) <= tol for t in grid) && t_min < 0.0 < t_max
        push!(grid, 0.0)
        sort!(grid)
    end

    flags = [is_feasible(current_point .+ t .* direction) for t in grid]
    intervals = Tuple{Float64, Float64}[]
    current_start = nothing

    for k in 1:(length(grid) - 1)
        t_left = grid[k]
        t_right = grid[k + 1]
        left_flag = flags[k]
        right_flag = flags[k + 1]

        if left_flag && isnothing(current_start)
            current_start = t_left
        end

        if left_flag != right_flag
            boundary = _locate_feasibility_transition(current_point, direction,
                t_left, left_flag, t_right, right_flag, is_feasible; tol, max_iter)

            if left_flag
                push!(intervals, (current_start, boundary))
                current_start = nothing
            else
                current_start = boundary
            end
        elseif left_flag && right_flag && k == length(grid) - 1
            push!(intervals, (current_start, t_right))
            current_start = nothing
        end
    end

    if !isempty(intervals)
        merged = Tuple{Float64, Float64}[intervals[1]]
        for interval in intervals[2:end]
            prev_lo, prev_hi = merged[end]
            lo, hi = interval
            if lo <= prev_hi + tol
                merged[end] = (prev_lo, max(prev_hi, hi))
            else
                push!(merged, interval)
            end
        end
        return merged
    end

    return intervals
end

"""
    _sample_from_interval_union(intervals; tol)

Samples uniformly from the one-dimensional measure of a union of intervals.
"""
function _sample_from_interval_union(intervals::Vector{Tuple{Float64, Float64}};
        tol::Float64=1e-8)
    lengths = [max(0.0, hi - lo) for (lo, hi) in intervals]
    total_length = sum(lengths)
    total_length <= tol && return nothing

    target = rand() * total_length
    cumulative = 0.0
    for (interval, length) in zip(intervals, lengths)
        cumulative += length
        if target <= cumulative
            lo, hi = interval
            return lo + rand() * (hi - lo)
        end
    end

    lo, hi = intervals[end]
    return lo + rand() * (hi - lo)
end

"""
    next_sample(current_point, box, is_feasible)

Generates the next sample using a full-line Hit-and-Run update. Along a random
ambient direction, the algorithm reconstructs the resolved feasible intervals of
the line and samples uniformly from their union. If those line intervals were
computed exactly, the uniform distribution over the full feasible set would be
stationary even for disconnected sets. In this implementation the interval
reconstruction is approximate, with resolution controlled by `max_step`.
"""
function next_sample(current_point::Vector{Float64}, box::Vector{Tuple{Float64, Float64}},
        is_feasible::Function; max_step::Float64=0.02, tol::Float64=1e-8,
        max_iter::Int=60, direction_basis::Union{Nothing, Matrix{Float64}}=nothing)
    dim = length(current_point)
    free_coordinates = findall(i -> (box[i][2] - box[i][1]) > tol, eachindex(box))
    isempty(free_coordinates) && return copy(current_point), 0

    while true
        direction = _random_direction(dim, direction_basis)
        if length(free_coordinates) < dim
            for i in eachindex(box)
                if !(i in free_coordinates)
                    direction[i] = 0.0
                end
            end
            dir_norm = norm(direction)
            dir_norm < tol && continue
            direction ./= dir_norm
        end
        t_min_box, t_max_box, success = _line_box_parameter_range(direction, current_point, box)
        if !success
            continue
        end

        intervals = _feasible_intervals_along_line(current_point, direction,
            t_min_box, t_max_box, is_feasible; max_step, tol, max_iter)
        isempty(intervals) && continue

        t_new = _sample_from_interval_union(intervals; tol)
        isnothing(t_new) && continue

        new_point = current_point .+ t_new .* direction
        if !is_feasible(new_point)
            continue
        end
        return new_point, 0
    end
end

export next_sample

end
