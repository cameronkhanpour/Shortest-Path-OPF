function _feasible_grid(bounds::Vector{Tuple{Float64, Float64}},
        is_feasible::Function; grid_points::Int=250)
    x_range = range(bounds[1][1], bounds[1][2], length=grid_points)
    y_range = range(bounds[2][1], bounds[2][2], length=grid_points)
    z = [Float64(is_feasible([x, y])) for y in y_range, x in x_range]
    return x_range, y_range, z
end

function _build_feasible_space_plot(bounds::Vector{Tuple{Float64, Float64}},
        is_feasible::Function,
        title_str::String;
        grid_points::Int=250)
    x_range, y_range, z = _feasible_grid(bounds, is_feasible; grid_points)

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
        xlabel=L"P_{G,2}\ [\mathrm{p.u.}]",
        ylabel=L"P_{G,3}\ [\mathrm{p.u.}]",
        legend=:topleft)
    plt.contour!(fig, x_range, y_range, z;
        levels=[0.5],
        color=boundary_color,
        linewidth=1.8,
        label="")

    return fig
end

function _save_plot_bundle(fig, base_name::String)
    png_name = "$(base_name).png"
    pdf_name = "$(base_name).pdf"
    plt.savefig(fig, png_name)
    plt.savefig(fig, pdf_name)
    println("Saved: $png_name")
    println("Saved: $pdf_name")
end

"""
    plot_feasible_space(sampled_points, bounds, is_feasible, case_name; metric, grid_points)

Plots the sampled feasible region for a 2D control space and saves both PNG and PDF copies.
"""
function plot_feasible_space(sampled_points::Vector{Vector{Float64}},
        bounds::Vector{Tuple{Float64, Float64}},
        is_feasible::Function,
        case_name::String;
        metric::Float64=NaN,
        grid_points::Int=250)
    if isempty(sampled_points)
        @warn "Cannot plot empty set of sampled points."
        return
    end

    dim = length(sampled_points[1])
    if dim != 2
        @warn "Feasible space plotting is only supported for 2D control spaces. Current: $(dim)D."
        return
    end

    println("Creating feasible space plot...")
    title_str = isnan(metric) ?
        "Feasible Space - $(case_name)" :
        "Feasible Space - $(case_name)   C(S) = $(round(metric, digits=4))"
    fig = _build_feasible_space_plot(bounds, is_feasible, title_str; grid_points)

    sample_color = get(CS.devon, 0.10)
    x_coords = [pt[1] for pt in sampled_points]
    y_coords = [pt[2] for pt in sampled_points]
    plt.scatter!(fig, x_coords, y_coords;
        label="Hit-and-Run samples",
        mc=sample_color,
        ms=4.5,
        ma=0.9,
        markerstrokewidth=0.25,
        markerstrokecolor=:white)

    date_str = Dates.format(Dates.today(), "yyyy-mm-dd")
    base_name = "$(date_str)_feasible_space_$(replace(case_name, "/" => "_"))"
    _save_plot_bundle(fig, base_name)
end

"""
    plot_shortest_path(p1, p2, all_samples, case_data, case_functions, case_name, plot_id,
                       is_feasible, bounds; grid_points)

Plots the shortest feasible path between two sampled points against the straight-line
reference, overlaid on the feasible space.
"""
function plot_shortest_path(p1::Vector{Float64}, p2::Vector{Float64},
        all_samples::Vector{Vector{Float64}},
        case_data,
        case_functions,
        case_name::String,
        plot_id::Int,
        is_feasible::Function,
        bounds::Vector{Tuple{Float64, Float64}};
        grid_points::Int=250)
    dim = length(p1)
    if dim != 2
        @warn "Path plotting is only supported for 2D control spaces."
        return
    end

    println("Computing shortest path for pair #$(plot_id)...")

    K = 19
    dt = 1 / (K + 1)
    tvec = collect(0.0:dt:1.0)
    tol_inner = 1e-6
    tol_pf = 1e-8
    iter_max_pf = 20

    v0, beta, _, _, path_data, _ = SP.get_feasible_path(case_data, p1, p2, tvec,
        1e-3, tol_inner, tol_pf, 100, iter_max_pf, 0.0125, 1e-6,
        case_data.qc_data.u0, false)

    if beta >= tol_inner
        @warn "Could not find a feasible path for pair #$(plot_id). Skipping plot."
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

    println("  Pair #$(plot_id): rho = $(round(rho, digits=4))")
    title_str = "Shortest Feasible Path - $(case_name)   Pair #$(plot_id), rho = $(round(rho, digits=3))"
    fig = _build_feasible_space_plot(bounds, is_feasible, title_str; grid_points)

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
    base_name = "$(date_str)_shortest_path_$(replace(case_name, "/" => "_"))_pair_$(plot_id)"
    _save_plot_bundle(fig, base_name)
end
