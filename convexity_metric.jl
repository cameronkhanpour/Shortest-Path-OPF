module ConvexityMetric

using Memento
using PowerModels
import PowerModels as PM
import Ipopt
Memento.setlevel!(getlogger(PowerModels), "error")

import Plots as plt
import ColorSchemes as CS
using LaTeXStrings
using Dates

include("power_flow_rect_qc.jl")
import .PowerFlowRectQC as QC
include("shortest_path_v2.jl")
import .ShortestPathOPF as SP

include("hit_and_run.jl")
import .HitAndRun
include("convexity_metric_plots.jl")

using LinearAlgebra
using Statistics
using Random

export sample_nonconvex_space, compute_ratio, compute_convexity_metric
export plot_feasible_space, plot_shortest_path, _create_feasibility_oracle, get_u_bounds
export analyze_case, analyze_cases


# --- Core Logic Functions ---

function _create_feasibility_oracle(case_file::String;
        ind_u::Union{Vector{Int64}, UnitRange{Int64}, Nothing}=nothing,
        u0_override=nothing)
    mpc = QC.load_case(case_file)
    case_qc = QC.compute_qc(mpc, false)
    if occursin("case9dongchan", case_file) && isnothing(ind_u)
        case_data = SP.load_case(case_file, case_qc, 2:3)
    else
        case_data = SP.load_case(case_file, case_qc, ind_u)
    end

    u0_ref = isnothing(u0_override) ? case_data.qc_data.u0 : u0_override
    tol_pf, iter_max_pf = 1e-8, 20
    case_functions = SP.make_functions(case_data, 1e-6, tol_pf, iter_max_pf, u0_ref)
    base_pf_seed = deepcopy(case_data.qc_data.x0)
    last_pf_seed = Ref(deepcopy(base_pf_seed))

    function is_feasible(u::Vector{Float64})
        for pf_seed in (last_pf_seed[], base_pf_seed, nothing)
            x, flag = case_data.pf_solver(u, u0_ref, tol_pf, iter_max_pf, pf_seed)
            flag == 0 && continue

            cons, pf_cons = case_functions.point_oracle(case_functions, u, x)
            if maximum(cons, init=-Inf) <= 1e-6 &&
                    maximum(abs.(pf_cons), init=0.0) <= tol_pf
                last_pf_seed[] = copy(x)
                return true
            end
        end
        return false
    end

    return is_feasible, case_data, case_functions, mpc, u0_ref
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

function _default_initial_point(is_feasible::Function, case_data, mpc)
    ind_u = case_data.ind_u
    u_start_full, _, u_end_full, _ = SP.get_endpoints(mpc, case_data, 1e-6,
        Random.MersenneTwister())
    u_start = u_start_full[ind_u] .+ case_data.PdQd[ind_u]
    u_end = u_end_full[ind_u] .+ case_data.PdQd[ind_u]
    u_mid = 0.5 .* (u_start .+ u_end)

    for candidate in (u_mid, u_end, u_start)
        if is_feasible(candidate)
            return candidate
        end
    end
    error("Could not find an initial feasible point for sampling.")
end

function _active_box_bound_count(point::Vector{Float64},
        bounds::Vector{Tuple{Float64, Float64}};
        tol::Float64=1e-5)
    count = 0
    for i in eachindex(point)
        if (bounds[i][2] - bounds[i][1]) <= tol
            continue
        end
        if (point[i] - bounds[i][1]) <= tol || (bounds[i][2] - point[i]) <= tol
            count += 1
        end
    end
    return count
end

function _max_pairwise_distance(points::Vector{Vector{Float64}})
    length(points) < 2 && return 0.0
    best = 0.0
    for i in 1:length(points), j in (i + 1):length(points)
        best = max(best, norm(points[i] - points[j]))
    end
    return best
end

function _project_to_box(point::Vector{Float64}, bounds::Vector{Tuple{Float64, Float64}})
    projected = copy(point)
    for i in eachindex(projected)
        projected[i] = clamp(projected[i], bounds[i][1], bounds[i][2])
    end
    return projected
end

function _find_sampling_seed(is_feasible::Function, case_data, mpc,
        bounds::Vector{Tuple{Float64, Float64}};
        opf_trials::Int=12)
    best_point = nothing
    best_active_count = typemax(Int)

    function consider_candidate(candidate)
        projected = _project_to_box(candidate, bounds)
        is_feasible(projected) || return
        active_count = _active_box_bound_count(projected, bounds)
        if active_count < best_active_count
            best_point = projected
            best_active_count = active_count
        end
        return
    end

    pilot_points = _randomized_opf_point_cloud(case_data, mpc; trials=opf_trials)
    for point in pilot_points
        consider_candidate(point)
    end
    for i in 1:length(pilot_points), j in (i + 1):length(pilot_points)
        consider_candidate(0.5 .* (pilot_points[i] .+ pilot_points[j]))
    end
    for i in 1:length(pilot_points), j in (i + 1):length(pilot_points), k in (j + 1):length(pilot_points)
        consider_candidate((pilot_points[i] .+ pilot_points[j] .+ pilot_points[k]) ./ 3.0)
    end

    if !isnothing(best_point)
        println("Using interiorised sampling seed with $(best_active_count) active box bounds.")
        return best_point
    end

    return _project_to_box(_default_initial_point(is_feasible, case_data, mpc), bounds)
end

function _opf_point_from_solution(case_data, opf_bus_solution)
    n = case_data.qc_data.n
    opf_x = zeros(Float64, 2 * n)
    for (bus, vir) in opf_bus_solution
        k = parse(Int, bus)
        opf_x[k] = vir["vr"]
        opf_x[n + k] = vir["vi"]
    end
    opf_u = (case_data.J_c + 0.5 .* sum(map((k) ->
        case_data.J_k[k] * opf_x[k], 1:(2 * n)))) * opf_x
    return opf_u[case_data.ind_u] .+ case_data.PdQd[case_data.ind_u]
end

function _randomized_opf_point_cloud(case_data, mpc;
        trials::Int=16,
        rng=Random.MersenneTwister(1))
    solver = PM.optimizer_with_attributes(Ipopt.Optimizer,
        "print_level" => 0,
        "constr_viol_tol" => 1e-6,
        "acceptable_constr_viol_tol" => 1e-6)

    points = Vector{Vector{Float64}}()
    for _ in 1:trials
        mpc_mod = deepcopy(mpc)
        for (_, gen) in mpc_mod["gen"]
            gen["cost"] = [0.0; randn(rng); 0.0]
        end

        jump_model = PM.instantiate_model(mpc_mod, PM.ACRPowerModel, PM.build_opf;
            ref_extensions=[])
        result = PM.optimize_model!(jump_model,
            relax_integrality=false,
            optimizer=solver,
            solution_processors=[])

        if !haskey(result, "solution") || !haskey(result["solution"], "bus")
            continue
        end
        push!(points, _opf_point_from_solution(case_data, result["solution"]["bus"]))
    end
    return points
end

function _pilot_direction_basis(points::Vector{Vector{Float64}};
        rel_singular_tol::Float64=1e-3)
    length(points) < 2 && return nothing

    mean_point = _sample_mean(points)
    centered = reduce(hcat, points) .- mean_point
    factor = svd(centered)
    isempty(factor.S) && return nothing

    scale = maximum(factor.S)
    scale <= 0.0 && return nothing
    keep = findall(>(rel_singular_tol * scale), factor.S)
    isempty(keep) && return nothing
    return Matrix(factor.U[:, keep])
end

# Core sampler: accepts pre-built oracle and case data to avoid redundant case loading.
function sample_nonconvex_space(is_feasible::Function, case_data, mpc, num_samples::Int;
        burn_in_samples::Int=500,
        thinning_factor::Int=20,
        initial_point=nothing,
        bounds=nothing,
        max_step::Float64=0.02,
        direction_basis::Union{Nothing, Matrix{Float64}}=nothing)
    ind_u = case_data.ind_u

    println("Finding initial feasible point...")
    bounding_box = isnothing(bounds) ? get_u_bounds(mpc, ind_u) : bounds
    if !isnothing(initial_point)
        current_point = _project_to_box(initial_point, bounding_box)
    elseif !isnothing(case_data.Q) && !isnothing(case_data.u_center)
        current_point = zeros(size(case_data.Q, 2))
    elseif occursin("case9dongchan", case_data.case_dir) && length(ind_u) == 2
        current_point = [0.5, 0.5]
    else
        current_point = _find_sampling_seed(is_feasible, case_data, mpc, bounding_box)
    end

    current_point = _project_to_box(current_point, bounding_box)

    if !is_feasible(current_point)
        error("Could not find an initial feasible point for sampling.")
    end
    if !isnothing(direction_basis)
        println("Using pilot direction basis of rank $(size(direction_basis, 2)).")
    end

    println("Starting burn-in period of $(burn_in_samples) samples...")
    for i in 1:burn_in_samples
        current_point, _ = HitAndRun.next_sample(current_point, bounding_box, is_feasible;
            max_step, direction_basis)
        current_point = _project_to_box(current_point, bounding_box)
        if i % 100 == 0
            println("  Burn-in: $(i)/$(burn_in_samples)")
        end
    end
    println("Burn-in finished.")

    sampled_points = [current_point]
    println("Collecting $(num_samples) samples (thinning: $(thinning_factor))...")
    for i in 1:(num_samples - 1)
        for _ in 1:thinning_factor
            current_point, _ = HitAndRun.next_sample(current_point, bounding_box, is_feasible;
                max_step, direction_basis)
            current_point = _project_to_box(current_point, bounding_box)
        end
        push!(sampled_points, current_point)
        if i % 10 == 0
            println("  Sample $(i+1)/$(num_samples)")
        end
    end
    println("Sampling finished. $(length(sampled_points)) points collected.")
    return sampled_points
end

# Convenience wrapper: loads the case internally (slower for repeated calls).
function sample_nonconvex_space(case_file::String, num_samples::Int;
        burn_in_samples::Int=500, thinning_factor::Int=20, max_step::Float64=0.02)
    is_feasible, case_data, _, mpc, _ = _create_feasibility_oracle(case_file)
    return sample_nonconvex_space(is_feasible, case_data, mpc, num_samples;
        burn_in_samples, thinning_factor, max_step)
end

"""
    compute_ratio(start_point, end_point, case_data, case_functions; K, u0_ref)

Computes the convexity ratio rho(x, y) = ||x - y|| / d_S(x, y) for a pair of
points using pre-built case infrastructure. Because the shortest path is solved
on a K-segment discretisation, this is an empirical approximation to the
continuous-space ratio.
"""
function compute_ratio(start_point::Vector{Float64}, end_point::Vector{Float64},
        case_data, case_functions; K::Int=19, u0_ref=case_data.qc_data.u0)
    try
        dim_u = length(start_point)
        n = case_data.qc_data.n
        dim_p = dim_u + 2 * n
        tol_inner = 1e-6
        tol_pf = 1e-8
        iter_max_pf = 20
        nu_0 = 1e-6
        mu = 1e-5

        dt = 1.0 / (K + 1)
        tvec = collect(0.0:dt:1.0)
        n_points = length(tvec)

        # Phase 1: find a feasible path.
        v0, beta, _, _, path_data, _ = SP.get_feasible_path(case_data, start_point, end_point,
            tvec, 1e-3, tol_inner, tol_pf, 100, iter_max_pf, 0.0125, nu_0, u0_ref, false)

        beta >= tol_inner && return NaN

        if isnothing(path_data)
            # The straight line is already feasible, so the ratio is exactly 1.
            return 1.0
        end

        # Shift gI values to match tol_inner, following barrier_case_time_v2.jl.
        if length(path_data) > 2
            for k in 1:(n_points - 2)
                path_data[1][k] .-= tol_inner
            end
        end

        # Phase 2: optimise to the shortest feasible path on the current discretisation.
        v, _, _, _, _ = SP.get_shortest_path(case_functions, tvec, v0;
            tol_outer=1e-3, iter_max=100, mu, nu_0, path_data)

        path_points = [v[((k-1)*dim_p+1):((k-1)*dim_p+dim_u)] for k in 1:n_points]
        pl = sum(norm(path_points[i+1] - path_points[i]) for i in 1:(n_points-1))
        pl0 = norm(end_point - start_point)

        (pl < 1e-12 || pl0 < 1e-12) && return NaN
        return clamp(pl0 / pl, 0.0, 1.0)
    catch e
        @warn "compute_ratio failed for a point pair: $e. Returning NaN."
        return NaN
    end
end

"""
    compute_convexity_metric(sampled_points, case_data, case_functions; K, u0_ref)

Computes the empirical convexity metric C(S) = mean_{i<j} rho(x_i, x_j). Any
pair for which the ratio cannot be computed is skipped and reported explicitly
so the returned average is honest about solver failures.
"""
function compute_convexity_metric(sampled_points::Vector{Vector{Float64}},
        case_data, case_functions; K::Int=19, u0_ref=case_data.qc_data.u0)
    num_points = length(sampled_points)
    num_points < 2 && return NaN

    total_pairs = div(num_points * (num_points - 1), 2)
    println("Computing convexity metric over $(total_pairs) pairs...")
    ratios = Float64[]
    skipped_pairs = 0

    pair_count = 0
    for i in 1:(num_points - 1), j in (i + 1):num_points
        pair_count += 1
        print("\r  Pair $(pair_count)/$(total_pairs)...")
        ratio = compute_ratio(sampled_points[i], sampled_points[j], case_data, case_functions;
            K, u0_ref)
        if isnan(ratio)
            skipped_pairs += 1
        else
            push!(ratios, ratio)
        end
    end
    println()

    isempty(ratios) && return NaN
    if skipped_pairs > 0
        @warn "Skipped $(skipped_pairs) of $(total_pairs) pairs while computing C(S). The reported value uses only successful pair evaluations."
    end
    println("Used $(length(ratios)) of $(total_pairs) pairs in the empirical average.")
    return mean(ratios)
end


# --- Slice Utilities ---

function _sample_mean(sampled_points::Vector{Vector{Float64}})
    sample_matrix = reduce(hcat, sampled_points)
    return vec(mean(sample_matrix, dims=2))
end

function _choose_slice_center(sampled_points::Vector{Vector{Float64}})
    mean_point = _sample_mean(sampled_points)
    center_idx = argmin(norm.(sampled_points .- Ref(mean_point)))
    return deepcopy(sampled_points[center_idx])
end

function _select_slice_positions(sampled_points::Vector{Vector{Float64}}; slice_dim::Int=2)
    dim = length(sampled_points[1])
    variances = [Statistics.var([pt[i] for pt in sampled_points]; corrected=false) for i in 1:dim]
    order = sortperm(variances, rev=true)
    return order[1:min(slice_dim, dim)], variances
end

function _slice_reference_u0(full_case_data, center_point::Vector{Float64})
    u0_ref = deepcopy(full_case_data.qc_data.u0)
    full_indices = full_case_data.ind_u
    u0_ref[full_indices] = center_point .- full_case_data.PdQd[full_indices]
    return u0_ref
end

function _slice_bounds_from_samples(sampled_points::Vector{Vector{Float64}},
        selected_positions::Vector{Int},
        natural_bounds::Vector{Tuple{Float64, Float64}};
        pad_fraction::Float64=0.15,
        min_pad::Float64=0.02)
    bounds = Vector{Tuple{Float64, Float64}}(undef, length(selected_positions))
    for (i, pos) in enumerate(selected_positions)
        values = [pt[pos] for pt in sampled_points]
        lo = minimum(values)
        hi = maximum(values)
        span = hi - lo
        pad = max(min_pad, pad_fraction * max(span, min_pad))
        nat_lo, nat_hi = natural_bounds[i]
        clipped_lo = max(nat_lo, lo - pad)
        clipped_hi = min(nat_hi, hi + pad)
        if clipped_lo >= clipped_hi
            bounds[i] = natural_bounds[i]
        else
            bounds[i] = (clipped_lo, clipped_hi)
        end
    end
    return bounds
end

function _control_label(mpc::Dict{String, Any}, u_idx::Int)
    n = length(mpc["bus"])
    if u_idx <= n
        return "P_G(bus $(u_idx)) [p.u.]"
    else
        return "V^2(bus $(u_idx - n)) [p.u.]"
    end
end

function _build_custom_2d_plot(bounds::Vector{Tuple{Float64, Float64}},
        is_feasible::Function,
        title_str::String,
        x_label::String,
        y_label::String;
        grid_points::Int=250)
    x_range = range(bounds[1][1], bounds[1][2], length=grid_points)
    y_range = range(bounds[2][1], bounds[2][2], length=grid_points)
    z = [Float64(is_feasible([x, y])) for y in y_range, x in x_range]

    plt.gr(size=(720, 540), dpi=180)
    region_color = get(CS.devon, 0.72)
    boundary_color = get(CS.devon, 0.98)

    fig = plt.heatmap(x_range, y_range, z;
        c=plt.cgrad([:white, region_color]),
        clims=(0.0, 1.0),
        colorbar=false,
        interpolate=false,
        aspect_ratio=:equal,
        framestyle=:box,
        grid=false,
        title=title_str,
        xlabel=x_label,
        ylabel=y_label,
        legend=:topleft)
    plt.contour!(fig, x_range, y_range, z;
        levels=[0.5],
        color=boundary_color,
        linewidth=1.8,
        label="")
    return fig
end

function _save_plot_bundle_local(fig, base_name::String)
    png_name = "$(base_name).png"
    pdf_name = "$(base_name).pdf"
    plt.savefig(fig, png_name)
    plt.savefig(fig, pdf_name)
    println("Saved: $png_name")
    println("Saved: $pdf_name")
end

function _plot_coordinate_slice_feasible_space(sampled_points::Vector{Vector{Float64}},
        bounds::Vector{Tuple{Float64, Float64}},
        is_feasible::Function,
        case_name::String,
        x_label::String,
        y_label::String;
        metric::Float64=NaN,
        grid_points::Int=250)
    title_str = isnan(metric) ?
        "Coordinate Slice - $(case_name)" :
        "Coordinate Slice - $(case_name)   C(S) = $(round(metric, digits=4))"
    fig = _build_custom_2d_plot(bounds, is_feasible, title_str, x_label, y_label; grid_points)

    sample_color = get(CS.devon, 0.10)
    plt.scatter!(fig, [pt[1] for pt in sampled_points], [pt[2] for pt in sampled_points];
        label="Hit-and-Run samples",
        mc=sample_color,
        ms=4.5,
        ma=0.9,
        markerstrokewidth=0.25,
        markerstrokecolor=:white)

    date_str = Dates.format(Dates.today(), "yyyy-mm-dd")
    base_name = "$(date_str)_coordinate_slice_$(replace(case_name, "/" => "_"))"
    _save_plot_bundle_local(fig, base_name)
end

function _plot_coordinate_slice_path(p1::Vector{Float64}, p2::Vector{Float64},
        all_samples::Vector{Vector{Float64}},
        case_data,
        case_functions,
        u0_ref,
        case_name::String,
        plot_id::Int,
        is_feasible::Function,
        bounds::Vector{Tuple{Float64, Float64}},
        x_label::String,
        y_label::String;
        grid_points::Int=250)
    println("Computing shortest path for coordinate-slice pair #$(plot_id)...")

    K = 19
    dt = 1 / (K + 1)
    tvec = collect(0.0:dt:1.0)
    tol_inner = 1e-6
    tol_pf = 1e-8
    iter_max_pf = 20

    v0, beta, _, _, path_data, _ = SP.get_feasible_path(case_data, p1, p2, tvec,
        1e-3, tol_inner, tol_pf, 100, iter_max_pf, 0.0125, 1e-6, u0_ref, false)

    if beta >= tol_inner
        @warn "Could not find a feasible path for coordinate-slice pair #$(plot_id). Skipping plot."
        return
    end

    n_points = length(tvec)
    local v
    if isnothing(path_data)
        v = v0
    else
        if length(path_data) > 2
            for k in 1:(n_points - 2)
                path_data[1][k] .-= tol_inner
            end
        end
        v, _, _, _, _ = SP.get_shortest_path(case_functions, tvec, v0; path_data=path_data)
    end

    dim_u = length(p1)
    n = case_data.qc_data.n
    dim_p = dim_u + 2 * n
    path_points = [v[((k-1)*dim_p+1):((k-1)*dim_p+dim_u)] for k in 1:length(tvec)]
    path_x = [pt[1] for pt in path_points]
    path_y = [pt[2] for pt in path_points]

    pl0 = norm(p2 - p1, 2)
    pl = sum(norm(path_points[i+1] - path_points[i], 2) for i in 1:(length(path_points)-1))
    rho = (pl > 1e-12) ? pl0 / pl : 1.0

    title_str = "Coordinate Slice Path - $(case_name)   Pair #$(plot_id), rho = $(round(rho, digits=3))"
    fig = _build_custom_2d_plot(bounds, is_feasible, title_str, x_label, y_label; grid_points)

    plt.scatter!(fig, [pt[1] for pt in all_samples], [pt[2] for pt in all_samples];
        label="",
        mc=:black,
        ms=2.5,
        ma=0.12,
        markerstrokewidth=0)
    plt.plot!(fig, [p1[1], p2[1]], [p1[2], p2[2]];
        label="Straight line",
        line=:dash,
        color=:gray45,
        lw=2.2)
    plt.plot!(fig, path_x, path_y;
        label="Shortest feasible path",
        color=get(CS.devon, 0.08),
        lw=3.0)
    plt.scatter!(fig, [p1[1]], [p1[2]];
        label="Start",
        mc=:goldenrod2,
        ms=7.5,
        markerstrokewidth=0.35,
        markerstrokecolor=:black)
    plt.scatter!(fig, [p2[1]], [p2[2]];
        label="End",
        mc=:tomato3,
        ms=7.5,
        markerstrokewidth=0.35,
        markerstrokecolor=:black)

    date_str = Dates.format(Dates.today(), "yyyy-mm-dd")
    base_name = "$(date_str)_coordinate_slice_path_$(replace(case_name, "/" => "_"))_pair_$(plot_id)"
    _save_plot_bundle_local(fig, base_name)
end

function _plot_random_pairs(sampled_points::Vector{Vector{Float64}},
        case_data,
        case_functions,
        u0_ref,
        case_name::String,
        is_feasible::Function,
        bounds::Vector{Tuple{Float64, Float64}},
        x_label::String,
        y_label::String;
        num_path_plots::Int=2,
        grid_points::Int=250)
    rng = Random.MersenneTwister(42)
    used_pairs = Set{Tuple{Int, Int}}()
    plots_made = 0
    while plots_made < num_path_plots &&
            length(used_pairs) < div(length(sampled_points) * (length(sampled_points) - 1), 2)
        i = rand(rng, 1:length(sampled_points))
        j = rand(rng, 1:length(sampled_points))
        i == j && continue
        pair = minmax(i, j)
        pair in used_pairs && continue
        push!(used_pairs, pair)
        plots_made += 1
        _plot_coordinate_slice_path(
            sampled_points[i], sampled_points[j],
            sampled_points, case_data, case_functions, u0_ref,
            case_name, plots_made, is_feasible, bounds, x_label, y_label;
            grid_points)
    end
end


# --- Analysis Workflow ---

function analyze_case(case_file::String;
        full_num_samples::Int=10,
        burn_in_samples::Int=250,
        thinning_factor::Int=10,
        max_step::Float64=0.02)
    case_name = splitext(basename(case_file))[1]
    println("\n=== $(case_name) ===")

    println("Loading full-space case context...")
    full_is_feasible, full_case_data, full_case_functions, mpc, full_u0_ref =
        _create_feasibility_oracle(case_file)
    full_dim = length(full_case_data.ind_u)
    println("Full control-space dimension: $(full_dim)")

    println("\n--- Full-Space Metric ---")
    full_samples = sample_nonconvex_space(
        full_is_feasible, full_case_data, mpc, full_num_samples;
        burn_in_samples, thinning_factor,
        max_step)
    all_samples_feasible = all(full_is_feasible.(full_samples))
    sample_spread = _max_pairwise_distance(full_samples)
    println("Sample diagnostics:")
    println("  all sampled points feasible: $(all_samples_feasible)")
    println("  max pairwise sample distance: $(round(sample_spread, digits=6))")

    full_metric = compute_convexity_metric(full_samples, full_case_data, full_case_functions;
        u0_ref=full_u0_ref)
    if !isnan(full_metric)
        println("Full-space C(S) = $(round(full_metric, digits=4))")
    end

    return (
        case_file=case_file,
        full_dim=full_dim,
        full_metric=full_metric,
        full_num_samples=full_num_samples,
        all_samples_feasible=all_samples_feasible,
        sample_spread=sample_spread,
    )
end

function analyze_cases(case_files::Vector{String}; kwargs...)
    results = NamedTuple[]
    for case_file in case_files
        push!(results, analyze_case(case_file; kwargs...))
    end
    return results
end

end  # module ConvexityMetric


# --- Main Execution ---
using .ConvexityMetric

function main(;
        case_files::Vector{String}=[
            "MATPOWER/pglib_opf_case60_c.m",
        ],
        full_num_samples::Int=10,
        burn_in_samples::Int=250,
        thinning_factor::Int=10,
        max_step::Float64=0.02)
    println("=== Convexity Metric Batch ===")
    println("Settings: full-space samples=$(full_num_samples), burn-in=$(burn_in_samples), thinning=$(thinning_factor), max_step=$(max_step)")

    results = ConvexityMetric.analyze_cases(case_files;
        full_num_samples,
        burn_in_samples,
        thinning_factor,
        max_step)

    println("\n=== Summary ===")
    for result in results
        println(splitext(basename(result.case_file))[1])
        println("  full-space dimension: $(result.full_dim)")
        println("  full-space C(S): $(isnan(result.full_metric) ? "NaN" : string(round(result.full_metric, digits=4)))")
        println("  all sampled points feasible: $(result.all_samples_feasible)")
        println("  max pairwise sample distance: $(round(result.sample_spread, digits=6))")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
