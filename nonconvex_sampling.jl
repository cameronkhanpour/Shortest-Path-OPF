module NonconvexSampling

using JuMP
using Ipopt
using SparseArrays
using LinearAlgebra
import PowerModelsAnnex as PMA

import MathOptInterface as MOI

include("power_flow_rect_qc.jl")
using .PowerFlowRectQC: parse_case, load_case   

# creates AC OPF JuMP model to give us the feasible space to sample from that particular case file 
function nonconvex_feasible_space(case_file::String)
    barrier_points = Vector{Vector{Float64}}()
    solutions = Vector{Vector{Float64}}()
    mpc = load_case(case_file)
    model = PMA.build_ac_opf(mpc)

    # gets all binding constraints of JuMP model
    binding_cons = ConstraintRef[]
    for (F, S) in list_of_constraint_types(model)  
        for con in all_constraints(model, F, S)
            if norm(dual(con)) > threshold || continue
                push!(binding_cons, con)
            end
        end
    end

end


# hit and run algorithm to generate uniform samples from AC OPF feasible space
function hit_and_run(is_inside::Bool, x_0::Float64, n_steps::Int64 = 10)

end


end