module ConvexityMetric

using Statistics, Random, Printf, LinearAlgebra

include("barrier_case_time_v2.jl")
include("power_flow_rect_qc.jl")
include("shortest_path_v2.jl")

import .PowerFlowRectQC: compute_qc
import .ShortestPathOPF: load_case, get_endpoints, run_case

export sample_nonconvex_space, compute_ratio, compute_convexity_metric

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

sample_points = sample_nonconvex_space(data, 5; ind_u = collect(2:3))
@printf("Sampled control points:\n")
for s in sample_points
    @printf("  %s\n", s)
end

using .ConvexityMetric: compute_convexity_metric
avg_ratio, ratios = compute_convexity_metric(sample_points; case_file = data)
println("Computed average convexity metric: ", avg_ratio)
