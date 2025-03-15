module NonconvexSampling

using JuMP
using Ipopt
using SparseArrays
using LinearAlgebra
import MathOptInterface as MOI

include("power_flow_rect_qc.jl")
using .PowerFlowRectQC: parse_case

# Barrier parameters (you can adjust these)
const μ_barrier = 1e-2   # weight on barrier term
const δ_barrier = 1e-4   # minimum allowed squared distance

"""
    build_opf_model(case_file::String, barrier_points::Vector{Vector{Float64}}, ind_u::Vector{Int}=nothing)

Builds a JuMP model whose decision variable is u.
If `ind_u` is provided then u is defined on those indices of qc_data.u0;
otherwise u is taken as the full qc_data.u0.
Here we use a dummy objective (sum(u.^2)) plus barrier terms.
(You can replace base_obj with a more meaningful objective if needed.)
"""
function build_opf_model(case_file::String, barrier_points::Vector{Vector{Float64}}, ind_u::Vector{Int}=nothing)
    # Parse the test case to obtain qc_data
    case, choosers, input, opf_u, opf_x, opf_fobj, mpc = parse_case(case_file)
    qc_data, _, _, _, _, _, _ = PowerFlowRectQC.compute_qc(case_file, false)
    if ind_u === nothing
        u_dim = length(qc_data.u0)
        indices = 1:u_dim
    else
        u_dim = length(ind_u)
        indices = ind_u
    end

    model = Model(Ipopt.Optimizer)
    set_silent(model)

    @variable(model, u[1:u_dim])
    # Set initial guess from qc_data.u0 restricted to indices
    for j in 1:u_dim
        set_start_value(u[j], qc_data.u0[indices][j])
    end

    # For this sampling example, we use a simple (dummy) base objective.
    # In practice you might use your ACOPF objective (transformed into the u space).
    @NLexpression(model, base_obj, sum(u[j]^2 for j=1:u_dim))
    @NLexpression(model, barrier_expr, sum( -log( sum((u[j] - bp[j])^2 for j=1:u_dim) - δ_barrier )
        for bp in barrier_points ) )
    @NLobjective(model, Min, base_obj + μ_barrier * barrier_expr)

    return model, qc_data
end

"""
    sample_nonconvex_space(case_file::String, num_samples::Int; ind_u::Vector{Int}=nothing)

Samples the nonconvex feasible space by repeatedly solving the u–space NLP.
For each sample the current solution is added as a barrier term.
"""
function sample_nonconvex_space(case_file::String, num_samples::Int; ind_u::Vector{Int}=nothing)
    barrier_points = Vector{Vector{Float64}}()
    solutions = Vector{Vector{Float64}}()
    
    for i in 1:num_samples
        println("Sampling iteration $i ...")
        model, qc_data = build_opf_model(case_file, barrier_points, ind_u)
        optimize!(model)
        status = termination_status(model)
        if status != MOI.OPTIMAL && status != MOI.LOCALLY_SOLVED
            println("Solver did not converge at iteration $i, status: $status")
            break
        end
        sol = [value(model[:u][j]) for j in 1:length(model[:u])]
        push!(solutions, sol)
        push!(barrier_points, sol)
        println("Iteration $i: Obtained a new sample; barrier updated.")
    end

    return solutions
end

# if abspath(PROGRAM_FILE) == @__FILE__
#     # Example: sample 5 points using the control variable space.
#     case_file = "MATPOWER/case9dongchan.m"
#     # For case9dongchan.m we want a 2D u: choose indices 2:3.
#     num_samples = 5
#     sols = sample_nonconvex_space(case_file, num_samples; ind_u = 2:3)
#     println("Collected $(length(sols)) solution samples.")
#     for (i, s) in enumerate(sols)
#         println("Sample $i: ", s)
#     end
# end

end  # module NonconvexSampling