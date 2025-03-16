module ConvexityMetric

using Statistics, Random, Printf, LinearAlgebra

# Include necessary repository files.
include("barrier_case_time_v2.jl")
include("power_flow_rect_qc.jl")
include("shortest_path_v2.jl")

# Import helper functions from the repository modules.
import .PowerFlowRectQC: compute_qc
import .ShortestPathOPF: load_case, get_endpoints, run_case

export sample_nonconvex_space, compute_ratio, compute_convexity_metric


function sample_nonconvex_space(case_file::String, num_samples::Int; ind_u::Vector{Int}=nothing)
    if case_file == "MATPOWER/case9dongchan.m" && ind_u === nothing
        ind_u = 2:3
    end

    case_qc = compute_qc(case_file, false)
    case_data = load_case(case_file, case_qc, ind_u)
    
    tol_inner = 1e-6

    u_start_full, x_start, u_end_full, x_end = get_endpoints(case_data.mpc, case_data, tol_inner)
    initial_sample = u_end_full[ind_u]
    samples = [initial_sample]
    @printf("Initial sample: %s\n", initial_sample)

    base_u = case_data.qc_data.u0[ind_u]
    
    for i in 2:num_samples
        repulsion = zeros(length(base_u))
        for p in samples
            diff = base_u .- p
            d = norm(diff)
            d = d < 1e-3 ? 1e-3 : d
            repulsion .+= diff ./ d
        end
        repulsion .*= 0.1 

        new_u_start = base_u .+ repulsion .+ 0.05 * randn(length(base_u))
        new_u_end   = base_u .- repulsion .+ 0.05 * randn(length(base_u))
        @printf("Iteration %d: Proposed endpoints:\n  new_u_start = %s\n  new_u_end   = %s\n", i, new_u_start, new_u_end)
        
        new_u_start_full = copy(case_data.qc_data.u0)
        new_u_end_full   = copy(case_data.qc_data.u0)
        new_u_start_full[ind_u] = new_u_start
        new_u_end_full[ind_u]   = new_u_end

        _ = run_case(case_file; u_start=new_u_start_full, u_end=new_u_end_full)
        
        u_start_full_new, x_start_new, u_end_full_new, x_end_new = get_endpoints(case_data.mpc, case_data, tol_inner, rng)
        new_sample = u_end_full_new[ind_u]
        @printf("Iteration %d: New sample: %s\n", i, new_sample)
        push!(samples, new_sample)
    end

    return samples
end

function compute_ratio(u_start::Vector{Float64}, u_end::Vector{Float64}, case_file::String)
    max_violation, max_violation_after, path_diff_pc, fobj_diff_pc, total_time =
         run_case(case_file; u_start=u_start, u_end=u_end)
    if isnan(fobj_diff_pc)
        println("Warning: run_case returned NaN for endpoints.")
        return NaN
    end
    return 1.0 + fobj_diff_pc / 100.0
end


function compute_convexity_metric(sample_points::Vector{Vector{Float64}}; case_file::String)
    n_samples = length(sample_points)
    if n_samples < 2
        error("At least two sample points are required to compute the convexity metric.")
    end
    ratios = Float64[]
    for i in 1:(n_samples-1)
        for j in (i+1):n_samples
            @printf("Computing ratio for pair (%d, %d)...\n", i, j)
            r = compute_ratio(sample_points[i], sample_points[j], case_file)
            if !isnan(r)
                push!(ratios, r)
                @printf("Pair (%d, %d): ratio = %.4f\n", i, j, r)
            else
                @printf("Pair (%d, %d): ratio is NaN, skipping.\n", i, j)
            end
        end
    end
    avg_ratio = mean(ratios)
    @printf("Average convexity metric (average ratio) = %.4f\n", avg_ratio)
    return avg_ratio, ratios
end

end  

# --- Run setup (example) ---

data = "MATPOWER/case9dongchan.m"
include("nonconvex_sampling.jl")
using .NonconvexSampling: sample_nonconvex_space

sample_points = sample_nonconvex_space(data, 5; ind_u = 2:3)
@printf("Sampled control points:\n")
for s in sample_points
    @printf("  %s\n", s)
end

using .ConvexityMetric: compute_convexity_metric
avg_ratio, ratios = compute_convexity_metric(sample_points; case_file = data)
println("Computed average convexity metric: ", avg_ratio)
