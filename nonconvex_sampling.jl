module NonconvexSampling

using JuMP
using Ipopt
using SparseArrays
using LinearAlgebra
import MathOptInterface as MOI

include("power_flow_rect_qc.jl")
using .PowerFlowRectQC: parse_case   


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


# test
case_file = "MATPOWER/case9dongchan.m"

num_samples = 5
sols = sample_nonconvex_space(case_file, num_samples; ind_u = collect(2:3))
println("Collected $(length(sols)) solution samples.")
for (i, s) in enumerate(sols)
    println("Sample $i: ", s)
end

end