using LinearAlgebra, Random, StatProfilerHTML, Profile, LaTeXStrings, Printf, JuMP 

import PowerModels as PM, Memento, Ipopt as Ipopt, Logging, SparseArrays

data = PM.parse_file("MATPOWER/case9mod.m")

PM.silence()

PM.print_summary(data)