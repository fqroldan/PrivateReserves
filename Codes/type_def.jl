using QuantEcon

mutable struct SOEres{K, Kd, Ktot}
	pars::Dict{Symbol, Float64}
	opt::Dict{Symbol, Bool}
	gr::Dict{Symbol, Vector{Float64}}
	prob::Dict{Symbol, Matrix{Float64}}

	v::Dict{Symbol, Array{Float64, K}}
	ϕ::Dict{Symbol, Array{Float64, Kd}}

	eq::Dict{Symbol, Array{Float64, Kd}}
	gov::Dict{Symbol, Array{Float64, Ktot}}
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
a 
	πLH=0.15,		# Transition prob for shock to spreads
	πHL=0.8,		# Transition prob for shock to spreads
	ψ=15,			# Inverse exposure of foreigners to domestic shock

	ϖ=0.6,			# Relative weight of nontradables
	η=1/0.83-1,		# Elasticity of substitution btw T and N

	wbar=0.95,		# Wage rigidity
	r=1.04^.25-1,	# Risk-free rate

	ρy=0.84,		# AR(1) for TFP in tradable sector
	σy=0.045,		# AR(1) for TFP in tradable sector

	α=0.75, 		# Curvature of production function
	Nb = 9,
	Na = 9,
	Nz = 5
	)
	
	ρy, σy = quarterlize_AR1(ρy, σy)
	μy = -0.5 * σy^2

	bgrid = range(-2,2,length=Nb)
	agrid = range(0,2,length=Na)

	zchain = tauchen(Nz, ρy, σy, 0, 1)
	zgrid, Pz = zchain.state_values, zchain.p

	Pϵ = [πLH 1-πLH; 1-πHL πHL]
	Nϵ = length(Pϵ[1,:])
	ϵgrid = range(0,1,length=Nϵ)

	R = zeros(Nb, Na, Nz, Nϵ)
	D = zeros(Nb, Na, Nz, Nϵ)
	V = zeros(Nb, Na, Nz, Nϵ)

	bp = ones(Nb, Na, Nz, Nϵ, 2) * 0.0
	ap = ones(Nb, Na, Nz, Nϵ, 2) * 0.5
	cc = ones(Nb, Na, Nz, Nϵ, 2) * 0.5
	
	cT = ones(Nb, Na, Nz, Nϵ, 2) * 0.5
	cN = ones(Nb, Na, Nz, Nϵ, 2) * 0.5
	hp = ones(Nb, Na, Nz, Nϵ, 2) * 0.5
	output = ones(Nb, Na, Nz, Nϵ, 2) * 0.5
	CA = ones(Nb, Na, Nz, Nϵ, 2) * 0.5

	repay = ones(Nb, Na, Nz, Nϵ, Nz, Nϵ)


	pars = Dict{Symbol, Float64}(:β=>β, :γ=>γ, :κ=>κ, :δ=>δ, :πLH=>πLH, :πHL=>πHL, :ψ=>ψ, :ϖ=>ϖ, :η=>η, :wbar=>wbar, :r=>r, :ρy=>ρy, :σy=>σy, :α=>α)
	opt = Dict{Symbol, Bool}()
	gr = Dict{Symbol, Vector{Float64}}(:b=>bgrid, :a=>agrid, :z=>zgrid, :ϵ=>ϵgrid)
	prob = Dict(:z=>Pz, :ϵ=>Pϵ)

	v = Dict(:R=>R, :D=>D, :V=>V)
	ϕ = Dict(:a=>ap, :b=>bp, :c=>cc)

	eq = Dict(:cT=>cT, :cN=>cN, :labor=>hp, :output=>output, :CA=>CA)
	gov = Dict(:repay => repay)

	return SOEres{4,5,6}(pars, opt, gr, prob, v, ϕ, eq, gov)
end
