using QuantEcon

mutable struct SOEres{K, Kd, Ktot}
	pars::Dict{Symbol, Float64}
	opt::Dict{Symbol, Bool}
	gr::Dict{Symbol, Vector{Float64}}
	prob::Dict{Symbol, Matrix{Float64}}

	v::Dict{Symbol, Array{Float64, K}}
	ϕ::Dict{Symbol, Array{Float64, Kd}}

	eq::Dict{Symbol, Array{Float64, Kd}}
	# gov::Dict{Symbol, Array{Float64, Ktot}}
end

function quarterlize_AR1(ρ, σ)
	ρ4 = ρ^0.25
	σ4 = sqrt(  σ^2 / ( 1 + ρ4^2 + ρ4^4 + ρ4^6 )  )
	return ρ4, σ4
end

function SOEres(;
	β=1.11^-.25,	# Discount factor
	γ=2.273,		# Risk aversion

	κ=0.185,		# Scale parameter for Extreme Value default choice shock
	δ=0.2845,		# Decay rate of government bonds
	
	πLH=0.15,		# Transition prob for shock to spreads
	πHL=0.8,		# Transition prob for shock to spreads
	ψ=15,			# Inverse exposure of foreigners to domestic shock
	θ=.04167,		# Reentry probability
	ℏ=.4,			# Haircut on default
	Δ=.1,			# Productivity loss in default

	ϖ=0.6,			# Relative weight of nontradables
	η=1/0.83-1,		# Elasticity of substitution btw T and N

	wbar=0.95,		# Wage rigidity
	r=1.04^.25-1,	# Risk-free rate

	ρz=0.84,		# AR(1) for TFP in tradable sector
	σz=0.045,		# AR(1) for TFP in tradable sector

	α=0.75, 		# Curvature of production function
	Nb = 9,
	Na = 9,
	Nz = 5
	)
	
	ρz, σz = quarterlize_AR1(ρz, σz)
	μz = -0.5 * σz^2

	bgrid = range(0,0.8,length=Nb)
	agrid = range(0,0.8,length=Na)

	zchain = tauchen(Nz, ρz, σz, μz, 1)
	zgrid, Pz = zchain.state_values, zchain.p

	Pν = [πLH 1-πLH; 1-πHL πHL]
	Nν = length(Pν[1,:])
	νgrid = range(0,1,length=Nν)

	R = zeros(Nb, Na, Nz, Nν)
	D = zeros(Nb, Na, Nz, Nν)
	V = zeros(Nb, Na, Nz, Nν)
	def = zeros(Nb, Na, Nz, Nν)

	bp = ones(Nb, Na, Nz, Nν, 2) * 0.5
	ap = ones(Nb, Na, Nz, Nν, 2) * 0.5
	cc = ones(Nb, Na, Nz, Nν, 2) * 0.5
	
	cT = ones(Nb, Na, Nz, Nν, 2) * 0.5
	cN = ones(Nb, Na, Nz, Nν, 2) * 0.5
	hp = ones(Nb, Na, Nz, Nν, 2) * 0.5
	output = ones(Nb, Na, Nz, Nν, 2) * 0.5
	CA = ones(Nb, Na, Nz, Nν, 2) * 0.5
	qb = ones(Nb, Na, Nz, Nν, 2) * exp(-r)
	qa = ones(Nb, Na, Nz, Nν, 2) * exp(-r)

	repay = ones(Nb, Na, Nz, Nν, Nz, Nν)


	pars = Dict{Symbol, Float64}(:β=>β, :γ=>γ, :κ=>κ, :δ=>δ, :πLH=>πLH, :πHL=>πHL, :ψ=>ψ, :ϖT=>ϖ, :ϖN=>1-ϖ, :η=>η, :wbar=>wbar, :r=>r, :ρz=>ρz, :σz=>σz, :α=>α, :θ=>θ, :ℏ=>ℏ, :Δ=>Δ, :ρz=>ρz, :σz=>σz, :μz=>μz)
	opt = Dict{Symbol, Bool}()
	gr = Dict{Symbol, Vector{Float64}}(:b=>bgrid, :a=>agrid, :z=>zgrid, :ν=>νgrid, :def=>0:1)
	prob = Dict(:z=>Pz, :ν=>Pν)

	v = Dict(:R=>R, :D=>D, :V=>V, :def=>def)
	ϕ = Dict(:a=>ap, :b=>bp, :c=>cT)

	eq = Dict(:cT=>cT, :cN=>cN, :labor=>hp, :output=>output, :CA=>CA, :qb=>qb, :qa=>qa)
	# gov = Dict(:repay => repay)

	return SOEres{4,5,6}(pars, opt, gr, prob, v, ϕ, eq)#, gov)
end

N(sr::SOEres, sym) = length(sr.gr[sym])

agg_grid(sr::SOEres) = gridmake(1:N(sr,:b), 1:N(sr,:a), 1:N(sr,:z), 1:N(sr,:ν))

statenames(sr::SOEres) = [:b,:a,:z,:ν]

S(sr::SOEres, jvec)  = Dict(sym => sr.gr[sym][jvec[jj]] for (jj, sym) in enumerate(statenames(sr)))
S_index(sr::SOEres, jvec) = Dict(sym => [jvec[jj]] for (jj, sym) in enumerate(statenames(sr)))
S_vec(sr::SOEres, s) = [s[key] for key in statenames(sr)]


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

function price_index(sr::SOEres, pN, pT = 1)
	ϖN, ϖT, η = [sr.pars[key] for key in [:ϖN, :ϖT, :η]]

	return (ϖN^(1/(1+η)) * pN^(η/(1+η)) + ϖT^(1/(1+η)) * pT^(η/(1+η)))^((1+η)/η)
end