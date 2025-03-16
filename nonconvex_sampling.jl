module NonconvexSampling

using JuMP
using Ipopt
using SparseArrays
using LinearAlgebra
import PowerModelsAnnex as PMA

import MathOptInterface as MOI

include("power_flow_rect_qc.jl")
using .PowerFlowRectQC: parse_case, load_case   

function sample_nonconvex_space(case_file::String)
    barrier_points = Vector{Vector{Float64}}()
    solutions = Vector{Vector{Float64}}()
    mpc = load_case(case_file)
    model = PMA.build_ac_opf(mpc)
    return model
end
end

# test
using .NonconvexSampling: sample_nonconvex_space
case_file = "MATPOWER/case9dongchan.m"
println(sample_nonconvex_space(case_file))