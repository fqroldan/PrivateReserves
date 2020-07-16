function iter_simul!(pp::Path{T}, sr::SOEres, t, itp_ϕb, itp_ϕa, itp_def, itp_qb, itp_qa, itp_pz) where T

	jζ = Int(getfrompath(pp, t, :jζ))
	jz = Int(getfrompath(pp, t, :jz))
	jν = Int(getfrompath(pp, t, :jν))

	def = def_state(sr, jζ)

	ζt = sr.gr[:def][jζ]
	zt = sr.gr[:z][jz]
	νt = sr.gr[:ν][jν]

	pz = itp_pz(zt, sr.gr[:z])
	pν = sr.prob[:ν][jν,:]

	at = getfrompath(pp, t, :a)
	bt = getfrompath(pp, t, :b)

	state = Dict(:a=>at, :b=>bt, :z=>zt, :ν=>νt)

	bp = itp_ϕb(bt, at, zt, νt, ζt)
	ap = itp_ϕa(bt, at, zt, νt, ζt)

	qav = itp_qa(bp,ap,zt,νt,ζt)
	eqm_t = get_eqm(sr, bp, ap, state, pz, pν, jζ, itp_def, qav, itp_qb)

	C = CES_aggregator(sr, eqm_t[:cT], eqm_t[:cN])

	pNt = price_nontradable(sr, eqm_t[:cT], eqm_t[:cN])

	qb = price_debt(sr, [bp,ap], zt, νt, pz, pν, itp_def, itp_qb, jdef=def)
	qa = exp(-sr.pars[:r])

	# Fill values of equilibrium at t
	fill_path!(pp, t, eqm_t)
	fill_path!(pp, t, Dict(:pN=>pNt, :C=>C, :qb=>qb, :qa=>qa))

	# Draw shocks for t+1
	probz = cumsum(pz)
	jzp = findfirst(probz .> rand())
	probν = cumsum(pν)
	jνp = findfirst(probν .> rand())

	zpv = sr.gr[:z][jzp]
	νpv = sr.gr[:ν][jνp]
	
	if def
		def_prob = (1-sr.pars[:θ])
	else
		def_prob = itp_def(bp,ap,zpv,νpv)
	end

	def_prime = rand() < def_prob

	if def_prime
		jζp = 1
	else
		jζp = 2
	end

	newdef = 0
	if def_prime && !def
		bp = (1-sr.pars[:ℏ]) * bp
		newdef = 1
	end
	reentry = 0
	if def && !def_prime
		reentry = 1
	end

	if t < T
		fill_path!(pp, t+1, Dict(:b=>bp, :a=>ap, :jζ=>jζp, :jz=>jzp, :jν=>jνp, :ζ=>sr.gr[:def][jζp], :newdef => newdef, :reentry => reentry, :z=>zpv, :ν=>νpv))
	end
	nothing
end

function simul(sr::SOEres, simul_length=4*10000, burn_in=4*1000)
	Random.seed!(1)
	T = simul_length + burn_in

	pp = Path(T=T)

	itp_ϕb = make_itp(sr, sr.ϕ[:b])
	itp_ϕa = make_itp(sr, sr.ϕ[:a])
	itp_def = make_itp(sr, sr.v[:def]);
	itp_qa  = make_itp(sr, sr.eq[:qa]);
	itp_qb  = make_itp(sr, sr.eq[:qb]);
	itp_pz = interpolate((sr.gr[:z], sr.gr[:z]), sr.prob[:z], Gridded(Linear()));

	fill_path!(pp, 1, Dict(:jζ=>2, :ζ=>sr.gr[:def][2], :jz=>1, :jν=>1, :a=>mean(sr.gr[:a]), :b=>mean(sr.gr[:b])))

	for jt in 1:T
		iter_simul!(pp, sr, jt, itp_ϕb, itp_ϕa, itp_def, itp_qb, itp_qa, itp_pz)
	end
	return trim_path(pp, burn_in)
end