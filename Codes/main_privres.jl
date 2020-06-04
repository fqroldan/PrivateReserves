println("Reserve Accumulation with a Private IIP")
using QuantEcon, Interpolations, Optim, Dates, Printf, Random, JLD, PlotlyJS, ColorSchemes, ORCA, JSON

include("type_def.jl")
include("handle_itps.jl")
include("comp_eqm.jl")
include("planner.jl")
include("reporting_routines.jl")

println("Constructor sr = SOEres()")
println("Main loop vfi!(sr; tol::Float64=1e-4, maxiter::Int64=500, verbose::Bool=false)")