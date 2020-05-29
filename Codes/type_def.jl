mutable struct SOEres
	pars::Dict{Symbol, Float64}
	opt::Dict{Symbol, Bool}
	gr::Dict{Symbol, Vector{Float64}}
	prob::Dict{Symbol, Matrix{Float64}}

	ϕ::Dict{Symbol, Array{Float64, Ktot}}
	v::Dict{Symbol, Array{Float64, Ktot}}

	eq::Dict{Symbol, Vector{Float64}}
	gov::Dict{Symbol, Vector{Float64}}
	LoM::Dict{Symbol, Array{Vector{Float64}, Kshocks}}
end
