print_save("Reserve Accumulation with a Private IIP")
using QuantEcon, Interpolations, LinearAlgebra, Optim, Dates, Printf, Random, JLD, PlotlyJS, ColorSchemes, JSON, DataFrames, GLM, NLopt, ForwardDiff

include("type_def.jl")
include("handle_itps.jl")
include("comp_eqm.jl")
include("planner.jl")
include("reporting_routines.jl")
include("plotting_routines.jl")
include("simul.jl")

print_save("Constructor sr = SOEres(; Euler = false, …)")
print_save("Main loop vfi!(sr; tol::Float64=1e-4, maxiter::Int64=500, verbose::Bool=false)")