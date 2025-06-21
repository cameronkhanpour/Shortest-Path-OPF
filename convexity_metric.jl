module ConvexityMetric

using Memento
using PowerModels 
Memento.setlevel!(getlogger(PowerModels), "error")

using Plots
using ColorSchemes

include("power_flow_rect_qc.jl")
import .PowerFlowRectQC as QC
include("shortest_path_v2.jl")
import .ShortestPathOPF as SP
include("barrier_case_time_v2.jl") 

include("hit_and_run.jl")
import .HitAndRun

using LinearAlgebra
using Statistics
using Random

export sample_nonconvex_space, compute_ratio, compute_convexity_metric
export plot_feasible_space, plot_shortest_path

# --- Plotting Functions ---

function plot_feasible_space(sampled_points::Vector{Vector{Float64}}, bounds::Vector{Tuple{Float64, Float64}}, is_feasible::Function, case_name::String; grid_points=100)
    if isempty(sampled_points)
        @warn "Cannot plot empty set of sampled points."
        return
    end

    dim = length(sampled_points[1])
    if dim != 2
        @warn "Feasible space plotting is only supported for 2D control spaces. The current space is $(dim)D."
        return
    end

    println("Creating a contour plot of the feasible space...")
    x_range = range(bounds[1][1], bounds[1][2], length=grid_points)
    y_range = range(bounds[2][1], bounds[2][2], length=grid_points)
    
    # Evaluate feasibility on the grid
    z = [is_feasible([x, y]) for y in y_range, x in x_range]
    
    # Create the plot
    p = contourf(x_range, y_range, z,
        levels=1, color=cgrad([:lightgray, :lightblue]),
        legend=false, aspect_ratio=:equal,
        title="Feasible Space for $case_name",
        xlabel="P_G2 (p.u.)",
        ylabel="P_G3 (p.u.)",
        framestyle=:box)

    # Plot the Hit-and-Run samples on top
    x_coords = [p[1] for p in sampled_points]
    y_coords = [p[2] for p in sampled_points]
    scatter!(p, x_coords, y_coords,
        label="Hit-and-Run Samples",
        mc=:navy,
        ms=3.5,
        ma=0.8,
        legend=:topleft,
        markerstrokewidth=0)
    
    filename = "feasible_space_$(replace(case_name, "/" => "_")).png"
    savefig(p, filename)
    println("Saved feasible space plot to: $filename")
end


function plot_shortest_path(p1::Vector{Float64}, p2::Vector{Float64}, all_samples::Vector{Vector{Float64}}, case_data, case_name::String, plot_id::Int, is_feasible::Function, bounds::Vector{Tuple{Float64, Float64}}; grid_points=100)
    dim = length(p1)
    if dim != 2
        @warn "Path plotting is only supported for 2D control spaces."
        return
    end

    println("Plotting shortest path for pair #$(plot_id) with feasible space background...")
    
    # 1. Create the base plot with the feasible space contour
    x_range = range(bounds[1][1], bounds[1][2], length=grid_points)
    y_range = range(bounds[2][1], bounds[2][2], length=grid_points)
    z = [is_feasible([x, y]) for y in y_range, x in x_range]
    
    p = contourf(x_range, y_range, z,
        levels=1, color=cgrad([:lightgray, :lightblue]),
        legend=false, aspect_ratio=:equal,
        title="Shortest Feasible Path for $case_name (Pair #$plot_id)",
        xlabel="P_G2 (p.u.)",
        ylabel="P_G3 (p.u.)",
        framestyle=:box)

    # 2. Add all other sampled points for context
    scatter!(p, [p[1] for p in all_samples], [p[2] for p in all_samples], 
        label="Other Samples", mc=:black, ms=2, ma=0.1, legend=:topleft, markerstrokewidth=0)

    # 3. Calculate the shortest path
    K = 19
    dt = 1 / (K + 1)
    tvec = collect(0.0:dt:1.0)
    tol_inner = 1e-6
    u0 = case_data.qc_data.u0
    tol_pf = 1e-8
    iter_max_pf = 20
    case_functions = SP.make_functions(case_data, tol_inner, tol_pf, iter_max_pf, u0)
    
    v0, beta, _, _, path_data, _ = SP.get_feasible_path(case_data, p1, p2, tvec, 1e-3, tol_inner, tol_pf, 100, iter_max_pf, 0.0125, 1e-6, u0, false)
    
    if beta >= tol_inner
        @warn "Could not find a feasible path for pair #$(plot_id). Skipping plot."
        return
    end

    v, _, _, _, _ = SP.get_shortest_path(case_functions, tvec, v0; path_data=path_data)
    
    dim_u = length(p1)
    dim_p = dim_u + 2 * case_data.qc_data.n
    path_points = [v[((k-1)*dim_p+1):((k-1)*dim_p+dim_u)] for k in 1:length(tvec)]
    path_x = [pt[1] for pt in path_points]
    path_y = [pt[2] for pt in path_points]

    # 4. Plot the paths and endpoints on top
    plot!(p, [p1[1], p2[1]], [p1[2], p2[2]], label="Straight Line", line=:dash, color=:gray, lw=2.5)
    plot!(p, path_x, path_y, label="Shortest Feasible Path", color=:crimson, lw=3)
    scatter!(p, [p1[1], p2[1]], [p1[2], p2[2]],
        label="Endpoints", mc=[:springgreen :darkorange], 
        ms=7, marker=:star5, markerstrokewidth=0.5, markerstrokecolor=:black)

    filename = "shortest_path_$(replace(case_name, "/" => "_"))_pair_$(plot_id).png"
    savefig(p, filename)
    println("Saved shortest path plot to: $filename")
end


# --- Core Logic Functions ---

function _create_feasibility_oracle(case_file::String)
    mpc = QC.load_case(case_file)
    case_qc = QC.compute_qc(mpc, false)
    if occursin("case9dongchan", case_file)
        case_data = SP.load_case(case_file, case_qc, 2:3)
    else
        case_data = SP.load_case(case_file, case_qc)
    end
    u0, tol_pf, iter_max_pf = case_data.qc_data.u0, 1e-8, 20
    case_functions = SP.make_functions(case_data, 1e-6, tol_pf, iter_max_pf, u0)
    
    function is_feasible(u::Vector{Float64})
        x, flag = case_data.pf_solver(u, u0, tol_pf, iter_max_pf)
        flag == 0 && return false
        cons, pf_cons = case_functions.point_oracle(case_functions, u, x)
        return maximum(cons, init=-Inf) <= 1e-6 && maximum(abs.(pf_cons), init=0.0) <= tol_pf
    end
    return is_feasible, case_data, mpc
end

function get_u_bounds(mpc::Dict{String, Any}, ind_u::Vector{Int64})
    n, dim_u = length(mpc["bus"]), length(ind_u)
    bounds = Vector{Tuple{Float64, Float64}}(undef, dim_u)
    bus_to_gen_map = Dict(gen["gen_bus"] => gen["index"] for (_, gen) in mpc["gen"])
    
    for (i, u_idx) in enumerate(ind_u)
        if u_idx <= n
            bus_idx = u_idx
            if haskey(bus_to_gen_map, bus_idx)
                gen = mpc["gen"]["$(bus_to_gen_map[bus_idx])"]
                bounds[i] = (gen["pmin"], gen["pmax"])
            else
                bounds[i] = (-Inf, Inf)
            end
        else
            bus = mpc["bus"]["$(u_idx - n)"]
            bounds[i] = (bus["vmin"]^2, bus["vmax"]^2)
        end
    end
    return bounds
end

function sample_nonconvex_space(case_file::String, num_samples::Int; mixing::Int=20)
    is_feasible, case_data, mpc = _create_feasibility_oracle(case_file)
    ind_u = case_data.ind_u
    
    println("Finding initial feasible point...")
    current_point = occursin("case9dongchan", case_file) ? [0.5, 0.5] :
        (SP.get_endpoints(mpc, case_data, 1e-6, Random.MersenneTwister())[1][ind_u] .+ case_data.PdQd[ind_u])

    if !is_feasible(current_point)
        error("Could not find an initial feasible point for sampling.")
    end

    bounding_box = get_u_bounds(mpc, ind_u)
    sampled_points = [current_point]
    
    println("Starting Hit-and-Run sampling for $(num_samples) samples (mixing factor: $(mixing))...")
    for i in 1:(num_samples - 1)
        for _ in 1:mixing
            current_point, _ = HitAndRun.next_sample(current_point, bounding_box, is_feasible)
        end
        push!(sampled_points, current_point)
        println("Collected sample $(i+1) / $(num_samples)")
    end
    println("Sampling finished.")
    return sampled_points
end

function compute_ratio(start_point::Vector{Float64}, end_point::Vector{Float64}, case_file::String)
    try
        _, _, _, fobj_diff_pc, _ = run_case(case_file; u_start=start_point, u_end=end_point, disable_msgs=true, K=19)
        isnan(fobj_diff_pc) ? NaN : (1.0 + fobj_diff_pc / 100.0)
    catch e
        @warn "An error occurred in run_case for a point pair: $e. Returning NaN."
        return NaN
    end
end

function compute_convexity_metric(sampled_points::Vector{Vector{Float64}}, case_file::String)
    num_points = length(sampled_points)
    num_points < 2 && return NaN
    
    total_pairs = Int(num_points * (num_points - 1) / 2)
    println("Computing convexity metric for $(total_pairs) pairs.")
    ratios = Float64[]
    
    pair_count = 0
    for i in 1:(num_points - 1), j in (i + 1):num_points
        pair_count += 1
        println("Processing pair $(pair_count) of $(total_pairs)...")
        ratio = compute_ratio(sampled_points[i], sampled_points[j], case_file)
        !isnan(ratio) && push!(ratios, ratio)
    end
    
    isempty(ratios) && return NaN
    return mean(ratios)
end

end 


# --- Main Test Execution ---
using .ConvexityMetric

function main()
    println("--- Starting Convexity Metric Test ---")
    case_file = "MATPOWER/case9dongchan.m"
    num_samples = 30
    num_path_plots = 3 

    # 1. Generate samples
    sampled_points = ConvexityMetric.sample_nonconvex_space(case_file, num_samples)
    println("\nGenerated $(length(sampled_points)) samples.")

    # 2. Create Plots
    println("\n--- Generating Plots ---")
    is_feasible, case_data, mpc = ConvexityMetric._create_feasibility_oracle(case_file)
    bounds = ConvexityMetric.get_u_bounds(mpc, case_data.ind_u)
    case_name = String(split(case_file, "/")[end])
    
    ConvexityMetric.plot_feasible_space(sampled_points, bounds, is_feasible, case_name)

    if length(sampled_points) >= 2
        for i in 1:min(num_path_plots, length(sampled_points) - 1)
            p1_idx, p2_idx = rand(1:length(sampled_points), 2)
            while p1_idx == p2_idx; p2_idx = rand(1:length(sampled_points)); end
            ConvexityMetric.plot_shortest_path(sampled_points[p1_idx], sampled_points[p2_idx], sampled_points, case_data, case_name, i, is_feasible, bounds)
        end
    end
    
    # 3. Compute the convexity metric
    metric = ConvexityMetric.compute_convexity_metric(sampled_points, case_file)
    
    println("\n--- Test Finished ---")
    if !isnan(metric)
        println("Calculated Convexity Metric for $(case_file): $(round(metric, digits=4))")
    else
        println("Failed to compute a convexity metric.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
