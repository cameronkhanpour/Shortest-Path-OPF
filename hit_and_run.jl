module HitAndRun

using LinearAlgebra, Random, Printf

"""
    _random_direction(dim::Int)

Generates a random unit vector in `dim` dimensions by sampling from a
multivariate normal distribution.
"""
function _random_direction(dim::Int)
    direction = randn(dim)
    return direction / norm(direction)
end

"""
    _line_box_intersection(direction, point, box)

Finds the intersection points of a line defined by `point + t*direction` with a
hyper-rectangle `box`. 
"""
#TODO maybe changing this to convex relaxation would be better  

function _line_box_intersection(direction::Vector{Float64}, point::Vector{Float64}, box::Vector{Tuple{Float64, Float64}})
    t_vals = Float64[]
    for i in 1:length(direction)
        if abs(direction[i]) > 1e-9
            t1 = (box[i][1] - point[i]) / direction[i]
            t2 = (box[i][2] - point[i]) / direction[i]
            push!(t_vals, t1)
            push!(t_vals, t2)
        end
    end

    valid_t = Float64[]
    for t in t_vals
        candidate_point = point .+ t .* direction
        is_inside = true
        for i in 1:length(direction)
            if !(box[i][1] - 1e-9 <= candidate_point[i] <= box[i][2] + 1e-9)
                is_inside = false
                break
            end
        end
        if is_inside
            push!(valid_t, t)
        end
    end

    if length(unique(round.(valid_t, digits=8))) < 2
        return nothing, nothing, false
    end
    
    t_min = minimum(valid_t)
    t_max = maximum(valid_t)
    p1 = point .+ t_min .* direction
    p2 = point .+ t_max .* direction
    
    return p1, p2, true
end

"""
    _find_boundary_by_stepping(current_point, direction, t_end, is_feasible; num_steps=200)

Finds the boundary of the feasible set by taking small, discrete steps along a direction
until an infeasible point is found.

Note: This method can be less computationally efficient than a bisection search,
especially if the feasible region is large along the search direction, as it
requires many calls to the expensive `is_feasible` function. The precision of
the boundary is also limited by the step size.

# Arguments
- `current_point`: The starting feasible point (where t=0).
- `direction`: The unit vector direction to search along.
- `t_end`: The maximum parameter value to search towards (from box intersection).
- `is_feasible`: The function that checks if a point is in the feasible set S.
- `num_steps`: The number of steps to take to cover the search interval.

# Returns
- The parameter `t` of the furthest known feasible point along the direction.
"""
function _find_boundary_by_stepping(current_point::Vector{Float64}, direction::Vector{Float64}, t_end::Float64, is_feasible::Function; num_steps=200)
    # Assumes the starting point (at t=0) is feasible.
    if num_steps <= 0
        error("num_steps must be positive.")
    end

    # The search starts from t=0 (current_point), which is known to be feasible.
    last_feasible_t = 0.0
    step_size = t_end / num_steps

    for i in 1:num_steps
        current_t = i * step_size
        
        # Check if the point at the current step is feasible
        if is_feasible(current_point + current_t * direction)
            # If it is, this is our new furthest known feasible point
            last_feasible_t = current_t
        else
            # If it's not feasible, we have crossed the boundary.
            # We stop and return the last point that was confirmed to be feasible.
            break
        end
    end

    return last_feasible_t
end


"""
    next_sample(current_point, box, is_feasible)

Generates the next sample using a line-search Hit-and-Run algorithm. This version
uses a stepping method to find the feasible boundaries.
"""
function next_sample(current_point::Vector{Float64}, box::Vector{Tuple{Float64, Float64}}, is_feasible::Function; num_steps=200)
    dim = length(current_point)
    
    while true 
        direction = _random_direction(dim)

        # 1. Find the intersection with the outer bounding box.
        p_box_1, p_box_2, success = _line_box_intersection(direction, current_point, box)
        if !success
            continue
        end
        
        # Convert intersection points to parameters 't' along the line.
        t_box_1 = dot(p_box_1 - current_point, direction)
        t_box_2 = dot(p_box_2 - current_point, direction)
        t_min_box = min(t_box_1, t_box_2)
        t_max_box = max(t_box_1, t_box_2)

        # 2. Use the stepping search to find the actual feasible boundaries.
        # Find positive boundary (t > 0)
        t_positive_bound = _find_boundary_by_stepping(current_point, direction, t_max_box, is_feasible; num_steps=num_steps)

        # Find negative boundary (t < 0). Search along the negative direction.
        # Note: -t_min_box is a positive value representing the distance to the negative boundary.
        t_negative_bound_in_pos_dir = _find_boundary_by_stepping(current_point, -direction, -t_min_box, is_feasible; num_steps=num_steps)
        t_negative_bound = -t_negative_bound_in_pos_dir

        # If the discovered feasible segment has negligible length, try a new direction.
        if abs(t_positive_bound - t_negative_bound) < 1e-9
            continue
        end

        # 3. Sample a new point uniformly from the found feasible line segment.
        lambda = rand()
        t_new = t_negative_bound + lambda * (t_positive_bound - t_negative_bound)
        
        new_point = current_point + t_new * direction

        return new_point, 0
    end
end

export next_sample

end
