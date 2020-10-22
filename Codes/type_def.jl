using QuantEcon, Distributions

mutable struct SOEres{K, Kd}
	pars::Dict{Symbol, Float64}
	opt::Dict{Symbol, Bool}
	gr::Dict{Symbol, Vector{Float64}}
	prob::Dict{Symbol, Matrix{Float64}}

	v::Dict{Symbol, Array{Float64, K}}
	ϕ::Dict{Symbol, Array{Float64, Kd}}

	eq::Dict{Symbol, Array{Float64, Kd}}
	# gov::Dict{Symbol, Array{Float64, Ktot}}
end

abstract type AbstractPath
end
mutable struct Path{T} <: AbstractPath
	data::Dict{Symbol, Vector{Float64}}
end
function Path(; T::Int64 = 1)
	data = Dict( key => Vector{Float64}(undef, T) for key in [:b, :a, :jζ, :jz, :jν, :ζ, :z, :ν, :pN, :C, :cT, :cN, :CA, :output, :labor, :qa, :qb, :spread, :newdef, :reentry, :π])
	return Path{T}(data)
end

periods(pv::Vector{T}) where T <: AbstractPath = sum([periods(pp) for pp in pv])
periods(p::Path{T}) where T = T

function check_periods(p::Path, t::Int64)
	0 < t <= periods(p) || throw("t out of bounds")
	nothing
end

getfrompath(p::Path, t::Int64, sym::Symbol) = p.data[sym][t]
getfrompath(p::Path, t::AbstractArray, sym::Symbol) = [p.data[sym][tv] for tv in t]
getfrompath(p::Path, t::Int) = Dict(key => p.data[key][t] for key in keys(p.data))
getfrompath(p::Path, sym::Symbol) = p.data[sym]
series(p::Path, sym::Symbol) = getfrompath(p,sym)
getmean(p::Path, sym::Symbol) = getmean([p], sym)
getmean(pv::Vector{T}, sym::Symbol) where T <: AbstractPath = mean(vcat([series(pp, sym) for pp in pv]...))

function fill_path!(p::Path, t::Int64, d::Dict=Dict())
	check_periods(p,t)
	missing_keys = 0
	for (key, val) in d
		if haskey(p.data, key)
			p.data[key][t] = val
		else
			missing_keys += 1
		end
	end
	
	if missing_keys > 0
		print_save("WARNING: $missing_keys missing keys")
	end
	nothing
end

function trim_path(p::Path{T}, t0::Int64) where T
	check_periods(p,t0)
	
	return Path{T-t0}(Dict(key => val[t0+1:end] for (key, val) in p.data))
end

function quarterlize_AR1(ρ, σ)
	ρ4 = ρ^0.25
	σ4 = sqrt(  σ^2 / ( 1 + ρ4^2 + ρ4^4 + ρ4^6 )  )
	return ρ4, σ4
end

function move_grids!(xgrid; xmin=0.0, xmax=1.0)
	xgrid[:] = xgrid[:] * (xmax-xmin) .+ xmin
 	nothing
 end

function SOEres(;
	β=1.06^-.25,	# Discount factor
	γ=2.273,		# Risk aversion

	σV=0.18,		# Scale parameter for Extreme Value default choice shock
	δ=0.2845,		# Decay rate of government bonds
	
	πLH=0.15,		# Transition prob for shock to spreads
	πHL=0.8,		# Transition prob for shock to spreads
	ψ=15,			# Inverse exposure of foreigners to domestic shock
	θ=.04167,		# Reentry probability
	ℏ=.4,			# Haircut on default
	Δ=.075,			# Productivity loss in default

	ϖ=0.55,			# Relative weight of nontradables
	η=1/0.83-1,		# Elasticity of substitution btw T and N

	wbar=0.7,		# Wage rigidity
	r=1.02^.25-1,	# Risk-free rate

	ρz=0.84,		# AR(1) for TFP in tradable sector
	σz=0.02,		# AR(1) for TFP in tradable sector

	α=0.75, 		# Curvature of production function
	Nb = 11,
	Na = 13,
	Nz = 9,
	bmax = 4.0,
	amax = 3.0
	)
	
	ρz, σz = quarterlize_AR1(ρz, σz)
	μz = -0.5 * σz^2

	bgrid = cdf.(Beta(2,1), range(0,1,length=Nb))
	agrid = cdf.(Beta(2,1), range(0,1,length=Na))
	move_grids!(bgrid, xmax=bmax)	
	move_grids!(agrid, xmax=amax)

	zchain = tauchen(Nz, ρz, σz, μz, 1)
	zgrid, Pz = zchain.state_values, zchain.p

	Pν = [πLH 1-πLH; 1-πHL πHL]
	Nν = length(Pν[1,:])
	νgrid = range(0,1,length=Nν)

	R = zeros(Nb, Na, Nz, Nν)
	D = zeros(Nb, Na, Nz, Nν)
	V = zeros(Nb, Na, Nz, Nν)
	def = zeros(Nb, Na, Nz, Nν)

	K = length(size(R))

	bp = ones(Nb, Na, Nz, Nν, 2) * mean(bgrid) * 0.75
	ap = ones(Nb, Na, Nz, Nν, 2) * mean(agrid)
	cc = ones(Nb, Na, Nz, Nν, 2) * 0.5
	
	cT = ones(Nb, Na, Nz, Nν, 2) * 0.5
	cN = ones(Nb, Na, Nz, Nν, 2) * 0.5
	hp = ones(Nb, Na, Nz, Nν, 2) * 0.5
	output = ones(Nb, Na, Nz, Nν, 2) * 0.5
	CA = ones(Nb, Na, Nz, Nν, 2) * 0.5
	qb = ones(Nb, Na, Nz, Nν, 2) * exp(-r)
	qa = ones(Nb, Na, Nz, Nν, 2) * exp(-r)

	repay = ones(Nb, Na, Nz, Nν, Nz, Nν)


	pars = Dict{Symbol, Float64}(:β=>β, :γ=>γ, :σV=>σV, :δ=>δ, :κ=>δ+exp(r)-1, :πLH=>πLH, :πHL=>πHL, :ψ=>ψ, :ϖT=>ϖ, :ϖN=>1-ϖ, :η=>η, :wbar=>wbar, :r=>r, :ρz=>ρz, :σz=>σz, :α=>α, :θ=>θ, :ℏ=>ℏ, :Δ=>Δ, :ρz=>ρz, :σz=>σz, :μz=>μz)
	opt = Dict{Symbol, Bool}()
	gr = Dict{Symbol, Vector{Float64}}(:b=>bgrid, :a=>agrid, :z=>zgrid, :ν=>νgrid, :def=>0:1)
	prob = Dict(:z=>Pz, :ν=>Pν)

	v = Dict(:R=>R, :D=>D, :V=>V, :def=>def)
	ϕ = Dict(:a=>ap, :b=>bp, :c=>cT, :cT=>cT)

	eq = Dict(:cT=>cT, :cN=>cN, :labor=>hp, :output=>output, :CA=>CA, :qb=>qb, :qa=>qa)
	# gov = Dict(:repay => repay)

	return SOEres{K,K+1}(pars, opt, gr, prob, v, ϕ, eq)#, gov)
end

N(sr::SOEres, sym) = length(sr.gr[sym])

agg_grid(sr::SOEres) = gridmake(1:N(sr,:b), 1:N(sr,:a), 1:N(sr,:z), 1:N(sr,:ν))

statenames(sr::SOEres) = [:b,:a,:z,:ν]

S(sr::SOEres, jvec)  = Dict(sym => sr.gr[sym][jvec[jj]] for (jj, sym) in enumerate(statenames(sr)))
S_index(sr::SOEres, jvec) = Dict(sym => [jvec[jj]] for (jj, sym) in enumerate(statenames(sr)))
S_vec(sr::SOEres, s) = [s[key] for key in statenames(sr)]

function corr(sr::SOEres, st)
	corr = ones(length(st))
	if any(statenames(sr).==:b)
		index_b = findfirst(statenames(sr).==:b)
		corr[index_b] *= (1-sr.pars[:ℏ])
		st_def = st .* corr
	else
		print("WARNING: NO VARIABLE CALLED :b IN TYPE")
		st_def = st
	end
	return st_def
end


def_state(sr::SOEres, jζ::Int64) = (jζ == 1)

function innov_z(sr::SOEres, zpv, zv)
	ρz, μz = [sr.pars[key] for key in [:ρz, :μz]]

	ϵpv = zpv - ρz * zv - (1-ρz) * μz
	return ϵpv
end

function CES_aggregator(sr::SOEres, cT, cN)
	ϖN, ϖT, η = [sr.pars[key] for key in [:ϖN, :ϖT, :η]]

	return (ϖN * cN^-η + ϖT * cT^-η)^(-1/η)
end

function price_nontradable(sr::SOEres, cT, cN, pT = 1)
	ϖN, ϖT, η = [sr.pars[key] for key in [:ϖN, :ϖT, :η]]

	pN = ϖN / (1-ϖT) * (cT/cN)^(1+η) * pT
end

function price_index(sr::SOEres, pN, pT = 1)
	ϖN, ϖT, η = [sr.pars[key] for key in [:ϖN, :ϖT, :η]]

	return (ϖN^(1/(1+η)) * pN^(η/(1+η)) + ϖT^(1/(1+η)) * pT^(η/(1+η)))^((1+η)/η)
end

function SDF(sr::SOEres, νv, ϵpv)
	ψ, σz, r = [sr.pars[sym] for sym in [:ψ, :σz, :r]]

	return exp(-r - νv * (ψ * ϵpv + 0.5 * ψ^2*σz^2))
end