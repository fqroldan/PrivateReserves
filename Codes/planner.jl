utility(sr::SOEres, c) = utility(c, sr.pars[:γ])
function utility(c, γ)
	cmin = 1e-8
	if c < cmin
		return utility(cmin,γ) + (c-cmin) * (cmin)^-γ
	else
		γ == 1 && return log(c)
		return c^(1-γ)/(1-γ)
	end
end

function prod_N(sr::SOEres, h)
	""" Computes production of nontradables at input h """
	yN = h^sr.pars[:α]
	return yN
end

function H(sr::SOEres, cT, w)
	""" Computes labor supply consistent with consumption of tradables + wage """
	α, η, ϖN, ϖT = [sr.pars[key] for key in [:α, :η, :ϖN, :ϖT]]

	return (ϖN/ϖT * α/w)^(1/(1+α*η)) * cT^((1+η)/(1+α*η))
end

function eq_h(sr::SOEres, cT)
	""" Computes labor supply consistent with consumption of tradables """
	Ls = 1

	h = H(sr, cT, sr.pars[:wbar])
	labor = min(h, Ls)
	return labor
end

bond_decay(sr::SOEres, jdef) = ifelse(jdef, 0.0, sr.pars[:δ])
output_T(sr::SOEres, state, jdef) = exp(state[:z]) * (1-sr.pars[:Δ]*jdef)

function price_debt(sr::SOEres, xp, zv, νv, pz, pν, itp_def, itp_q; jdef::Bool=false)
	""" Iterates once on the debt price using next period's state """
	δ, κC, ℏ, ψ, σz, r = [sr.pars[sym] for sym in [:δ, :κC, :ℏ, :ψ, :σz, :r]]
	qb = 0.0
	bp, ap = xp
	for (jzp, zpv) in enumerate(sr.gr[:z]), (jνp, νpv) in enumerate(sr.gr[:ν])
		prob = pz[jzp] * pν[jνp]

		ϵpv = innov_z(sr, zpv, zv)
		sdf = SDF(sr, νv, ϵpv)

		jζp = 1 # Default
		rep_default= (1-ℏ) * itp_q(bp*(1-ℏ),ap,zpv,νpv,sr.gr[:def][jζp])
		jζp = 2 # Repayment
		rep_normal = κC + (1-δ) * itp_q(bp,ap,zpv,νpv,sr.gr[:def][jζp])

		if jdef
			def_prob = sr.pars[:θ]
		else
			def_prob = itp_def(bp,ap,zpv,νpv)
		end

		qb += prob * sdf * (def_prob * rep_default + (1-def_prob) * rep_normal)
	end
	return qb
end

function budget_constraint_T(sr::SOEres, state, pz, pν, xp, itp_def, itp_q, qav, jdef::Bool)
	""" Computes consumption of T given state and choices of new debt and reserves """

	yT = output_T(sr, state, jdef)
	δv = bond_decay(sr, jdef)

	bv, av, zv, νv = [state[key] for key in [:b, :a, :z, :ν]]
	# qa = exp(-sr.pars[:r])

	bp, ap = xp

	debt_operations = 0.0
	if !jdef
		qb = price_debt(sr, xp, zv, νv, pz, pν, itp_def, itp_q)
		debt_operations = qb * (bp - (1-δv) * bv) - sr.pars[:κC] * bv
	end

	cT = yT + av - qav*ap + debt_operations
	return cT
end

function budget_constraint_agg(sr::SOEres, state, pz, pν, xp, itp_def, itp_q, qav, jdef::Bool)
	""" Computes aggregate of C at state and choices of new debt and reserves """
	cT = budget_constraint_T(sr, state, pz, pν, xp, itp_def, itp_q, qav, jdef)

	cT = max(0.0, cT)

	h = eq_h(sr, cT)
	yN = prod_N(sr, h)
	c = CES_aggregator(sr, cT, yN)

	return c
end

function value(sr::SOEres, state, pz, pν, bp, ap, itp_v, itp_vd, itp_def, itp_q, qav, jdef::Bool)
	""" Computes V given state and choices of consumption, debt, reserves """
	θ, ℏ, β = [sr.pars[sym] for sym in [:θ, :ℏ, :β]]

	# xp = Dict([(:b,bp), (:a,ap)])

	xp = [bp, ap]

	# bp = max(bp, minimum(sr.gr[:b]))
	# bp = min(bp, maximum(sr.gr[:b]))

	c = budget_constraint_agg(sr, state, pz, pν, xp, itp_def, itp_q, qav, jdef)
	ut = utility(sr, c)

	vp = 0.0
	for (jzp, zpv) in enumerate(sr.gr[:z]), (jνp, νpv) in enumerate(sr.gr[:ν])
		prob = pz[jzp] * pν[jνp]
		if jdef
			Vpv = θ * itp_v(xp..., zpv, νpv) + (1-θ) * itp_vd(xp..., zpv, νpv)
		else
			Vpv = itp_v(xp..., zpv, νpv)
		end

		vp += prob * Vpv
	end

	vt = ut + β * vp
	return vt
end

function opt_value_R(sr::SOEres, guess, state, pz, pν, itp_v, itp_vd, itp_def, itp_q, qav)
	""" Choose ap ≥ 0, bp ≥ 0 default with prob sr.gov[:repay] """
	jζ = 2 # IN REPAYMENT
	jdef = def_state(sr, jζ)
	jdef == false || throw(error("Wrong default state"))

	xguess = [guess[key][jζ] for key in [:b, :a]]
	xmin = [minimum(sr.gr[key]) for key in [:b, :a]]
	xmax = [maximum(sr.gr[key]) for key in [:b, :a]]
	xmax[1] = max(0.5*xmax[1], min(xmax[1], 2*state[:b]))

	obj_f(x) = -value(sr, state, pz, pν, x[1], x[2], itp_v, itp_vd, itp_def, itp_q, qav, jdef)
	res = Optim.optimize(obj_f, xmin, xmax, xguess, Fminbox(NelderMead()))

	if !Optim.converged(res)
		res = Optim.optimize(obj_f, xmin, xmax, xguess, Fminbox(GradientDescent()))
	end
	
	!Optim.converged(res) && println("WARNING: DIDN'T FIND SOL IN R")
	
	x_opt = res.minimizer
	bpv, apv = x_opt
	vR = -obj_f(x_opt)

	ϕ = Dict(:a=>apv, :b=>bpv)

	ccv = budget_constraint_agg(sr, state, pz, pν, x_opt, itp_def, itp_q, qav, jdef)
	ϕ[:c] = ccv
	return ϕ, vR
end

function opt_value_D(sr::SOEres, guess, state, pz, pν, itp_v, itp_vd, itp_def, itp_q, qav)
	""" Choose ap between 0 and a, bp = bv, reenter mkts w prob θ """
	bpv = state[:b]
	jζ = 1 # IN DEFAULT
	jdef = def_state(sr, jζ)
	jdef == true || throw(error("Wrong default state"))

	xguess = [guess[key][jζ] for key in [:a]]
	xmin = [minimum(sr.gr[key]) for key in [:a]]
	xmax = [maximum(sr.gr[key]) for key in [:a]]

	# obj_f(x) = -value(sr, state, pz, pν, bpv, x[1], itp_v, itp_vd, itp_def, itp_q, qav, jdef)
	# res = Optim.optimize(obj_f, xmin, xmax, xguess, Fminbox(GradientDescent()))
	obj_f(apv) = -value(sr, state, pz, pν, bpv, apv, itp_v, itp_vd, itp_def, itp_q, qav, jdef)
	res = Optim.optimize(obj_f, first(xmin), first(xmax), GoldenSection())

	!Optim.converged(res) && println("WARNING: DIDN'T FIND SOL IN D")

	x_opt = res.minimizer
	apv = first(x_opt)
	vD = -obj_f(x_opt)

	ϕ = Dict(:a=>apv, :b=>bpv)
	
	ccv = budget_constraint_agg(sr, state, pz, pν, [bpv, apv], itp_def, itp_q, qav, jdef)
	ϕ[:c] = ccv
	return ϕ, vD
end

function solve_optvalue(sr::SOEres, itp_v, itp_vd, itp_def, itp_q)
	""" Loop over states and solve """
	new_v = Dict(key => similar(val) for (key, val) in sr.v)
	new_ϕ = Dict(key => similar(val) for (key, val) in sr.ϕ)

	Jgrid = agg_grid(sr);
	Threads.@threads for js in 1:size(Jgrid,1)
		jv = Jgrid[js,:]

		jd = S_index(sr, jv)
		state = S(sr, jv)

		pz = sr.prob[:z][jd[:z],:]
		pν = sr.prob[:ν][jd[:ν],:]

		guess = Dict(key => val[jv...,:] for (key,val) in sr.ϕ)

		qaD, qaR = zeros(2)
		for jζ = 1:2
			if def_state(sr, jζ) == 1
				qaD = sr.eq[:qa][jv..., jζ]
			else
				qaR = sr.eq[:qa][jv..., jζ]
			end
		end

		""" New xp and values in repayment and default """
		ϕR, vR = opt_value_R(sr, guess, state, pz, pν, itp_v, itp_vd, itp_def, itp_q, qaR)
		ϕD, vD = opt_value_D(sr, guess, state, pz, pν, itp_v, itp_vd, itp_def, itp_q, qaD)

		""" Repayment probability as function of values """
		prob_rep = prob_extreme_value(sr, vR, vD)

		""" Save new values """
		new_v[:D][jv...] = vD
		new_v[:R][jv...] = vR
		for key in keys(sr.ϕ)
			for jζ in 1:2
				if def_state(sr, jζ)
					new_ϕ[key][jv...,jζ] = ϕD[key]
				else
					new_ϕ[key][jv...,jζ] = ϕR[key]
				end
			end
		end
	end

	return new_v, new_ϕ
end

function prob_extreme_value(sr::SOEres, vR, vD)
	""" Apply extreme-value shocks formula """
	κV = sr.pars[:κV]
	prob = exp(vD/κV) / (exp(vR/κV) + exp(vD/κV))
	if isnan(prob)
		if vR > vD
			return 0.0
		else
			return 1.0
		end
	end
	return prob
end

function update_def!(sr::SOEres)
	""" Computes default prob and value of entering period in repayment """
	""" Uses 'updated' sr.v[:D] and sr.v[:R] """
	itp_vd = make_itp(sr,sr.v[:D]);
	ℏ = sr.pars[:ℏ]

	Jgrid = agg_grid(sr);
	for js in 1:size(Jgrid,1)
		jv = Jgrid[js,:]

		state = S(sr, jv)
		st = [state[key] for key in statenames(sr)]
		st_def = corr(sr, st)

		vR = sr.v[:R][jv...]
		vD = itp_vd(st_def...)

		def_prob = prob_extreme_value(sr,vR,vD)

		sr.v[:def][jv...] = def_prob
		sr.v[:V][jv...] = def_prob * vD + (1-def_prob) * vR
	end
end

function vfi_iter(sr::SOEres)
	""" Interpolate values and prices to use as next period values """
	itp_v  = make_itp(sr, sr.v[:V]);
	itp_vd = make_itp(sr, sr.v[:D]);
	itp_def = make_itp(sr, sr.v[:def]);
	itp_q  = make_itp(sr, sr.eq[:qb]);

	new_v, new_ϕ = solve_optvalue(sr, itp_v, itp_vd, itp_def, itp_q);

	return new_v, new_ϕ
end

function update_sr!(y, new_y, upd_η = 1)
	for key in keys(y)
		y[key] = y[key] + upd_η * (new_y[key] - y[key])
	end
	nothing
end

function vfi!(sr::SOEres; tol::Float64=1e-4, maxiter::Int64=1000, verbose::Bool=false)
	""" Main Loop """
	iter, dist = 0, 1+tol
	avg_time = 0.0
	dist_v, dist_ϕ = zeros(2)

	upd_η = 0.25

	t0 = time()
	while dist > tol && iter < maxiter
		iter += 1

		old_q = copy(sr.eq[:qb]);
		""" Update debt prices (for use as next period prices) """
		update_q!(sr, verbose = verbose)
		dist_q = norm(sr.eq[:qb]-old_q) / max(1,norm(old_q))
		# sr.eq[:qb] = old_q + upd_ηq * (sr.eq[:qb] - old_q)
		norm_q = norm(sr.eq[:qb])

		t1 = time()
		""" Iterate on the value functions """
		new_v, new_ϕ = vfi_iter(sr);
		t = time() - t1

		avg_time = (avg_time * (iter - 1) + t) / iter

		dist_v = maximum([ norm(new_v[key] - sr.v[key]) / max(1,norm(sr.v[key])) for key in keys(sr.v) ])
		dist_ϕ = maximum([ norm(new_ϕ[key] - sr.ϕ[key]) / max(1,norm(sr.ϕ[key])) for key in keys(sr.ϕ) ])

		norm_v = maximum([norm(sr.v[key]) for key in keys(sr.v)])
		norm_ϕ = maximum([norm(sr.ϕ[key]) for key in keys(sr.ϕ)])
		
		dist = max(dist_v, dist_ϕ, dist_q)

		update_sr!(sr.v, new_v, upd_η)
		update_def!(sr) # To update [:def] and [:V] consistently with [:D], [:R]
		update_sr!(sr.ϕ, new_ϕ, upd_η)

		upd_η = max(upd_η * 0.995, 1e-2)

		# for key in keys(sr.v)
		# 	print("||$(key)|| = $(norm(sr.v[key]))\n")
		# end

		if verbose && iter % 1 == 0
			print("After $iter iterations (avg time = $(time_print(avg_time))), d(v,ϕ,q) = $(@sprintf("%0.3g",dist_v)), $(@sprintf("%0.3g",dist_ϕ)), $(@sprintf("%0.3g",dist_q)) \n")
			print("‖v,ϕ,q‖ = $(@sprintf("%0.3g",norm_v)), $(@sprintf("%0.3g",norm_ϕ)), $(@sprintf("%0.3g",norm_q)) \n")
			print("upd_η = $(@sprintf("%0.3g", upd_η))\n")
		end
	end

	update_eqm!(sr)

	if dist <= tol
		s = "Converged in $iter iterations (total time = $(time_print(time()-t0)))\n"
	else
		s = "After $iter iterations (total time = $(time_print(time()-t0))), d(v,ϕ) = $(@sprintf("%0.3g",dist_v)), $(@sprintf("%0.3g",dist_ϕ))\n"
	end
	verbose ? print(s) : nothing
	return dist
end
