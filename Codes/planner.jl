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
function uT(sr::SOEres, cT, yN)
	c = CES_aggregator(sr, cT, yN)

	return utility(sr,c)
end
uT_prime(sr, cT, yN) = ForwardDiff.derivative(x->uT(sr, x, yN), cT)

function prod_N(sr::SOEres, h, jdef)
	""" Computes production of nontradables at input h """

	yN = (1-sr.pars[:Δ]*jdef) * h^sr.pars[:α]
	# yN = h^sr.pars[:α]
	return yN
end

function H(sr::SOEres, cT, w, jdef)
	""" Computes labor supply consistent with consumption of tradables + wage """
	α, η, ϖN, ϖT = [sr.pars[key] for key in [:α, :η, :ϖN, :ϖT]]

	return (ϖN/ϖT * (1-sr.pars[:Δ]*jdef) * α/w)^(1/(1+α*η)) * cT^((1+η)/(1+α*η))
end

function eq_h(sr::SOEres, cT, jdef)
	""" Computes labor supply consistent with consumption of tradables """
	Ls = 1

	h = H(sr, cT, sr.pars[:wbar], jdef)
	labor = min(h, Ls)
	return labor
end

bond_decay(sr::SOEres, jdef) = ifelse(jdef, 0.0, sr.pars[:δ])
output_T(sr::SOEres, state, jdef) = exp(state[:z]) * (1-sr.pars[:Δ]*jdef)

function price_debt(sr::SOEres, xp, zv, νv, pz, pν, itp_def, itp_q; jdef::Bool=false)
	""" Iterates once on the debt price using next period's state """
	δ, κ, ℏ, ψ, σz, r = [sr.pars[sym] for sym in [:δ, :κ, :ℏ, :ψ, :σz, :r]]
	qb = 0.0
	bp, ap = xp
	for (jzp, zpv) in enumerate(sr.gr[:z]), (jνp, νpv) in enumerate(sr.gr[:ν])
		prob = pz[jzp] * pν[jνp]

		ϵpv = innov_z(sr, zpv, zv)
		sdf = SDF(sr, νv, ϵpv)

		jζp = 1 # Default
		rep_default= (1-ℏ) * itp_q(bp*(1-ℏ),ap,zpv,νpv,sr.gr[:def][jζp])
		jζp = 2 # Repayment
		rep_normal = κ + (1-δ) * itp_q(bp,ap,zpv,νpv,sr.gr[:def][jζp])

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
		debt_operations = qb * (bp - (1-δv) * bv) - sr.pars[:κ] * bv
	end

	cT = yT + av - qav*ap + debt_operations
	return cT
end

function budget_constraint_agg(sr::SOEres, state, pz, pν, xp, itp_def, itp_q, qav, jdef::Bool)
	""" Computes aggregate of C at state and choices of new debt and reserves """
	cT = budget_constraint_T(sr, state, pz, pν, xp, itp_def, itp_q, qav, jdef)

	cT = max(0.0, cT)

	h = eq_h(sr, cT, jdef)
	yN = prod_N(sr, h, jdef)
	c = CES_aggregator(sr, cT, yN)

	return c, cT
end

function value(sr::SOEres, state, pz, pν, bp, ap, itp_v, itp_vd, itp_def, itp_q, qav, jdef::Bool)
	""" Computes V given state and choices of consumption, debt, reserves """
	θ, ℏ, β = [sr.pars[sym] for sym in [:θ, :ℏ, :β]]

	# xp = Dict([(:b,bp), (:a,ap)])

	xp = [bp, ap]

	# bp = max(bp, minimum(sr.gr[:b]))
	# bp = min(bp, maximum(sr.gr[:b]))

	c, _ = budget_constraint_agg(sr, state, pz, pν, xp, itp_def, itp_q, qav, jdef)
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

function Euler_eq(sr::SOEres, state, pz, pν, bpv, apv, itp_ucT, itp_def, itp_q, qav, jdef::Bool)
	cT = budget_constraint_T(sr, state, pz, pν, [bpv, apv], itp_def, itp_q, qav, jdef)

	cT = max(0.0, cT)
	h = eq_h(sr, cT, jdef)
	yN = prod_N(sr, h, jdef)

	uT = uT_prime(sr, cT, yN)
	LHS = uT * qav

	EuT = 0.0
	for (jzp, zpv) in enumerate(sr.gr[:z]), (jνp, νpv) in enumerate(sr.gr[:ν])
		prob = pz[jzp] * pν[jνp]
		if jdef
			uTp = θ * itp_ucT(bpv, apv, zpv, νpv, 2) + (1-θ) * itp_ucT(bpv, apv, zpv, νpv, 1)
		else
			defprob = itp_def(bpv, apv, zpv, νpv)
			uTp = (1-defprob) * itp_ucT(bpv, apv, zpv, νpv, 2) + defprob * itp_ucT(bpv, apv, zpv, νpv, 1)
		end
		EuT += prob * uTp
	end

	RHS = β * EuT
	return RHS - LHS
end

function opt_value_R(sr::SOEres, guess, state, pz, pν, itp_ucT, itp_v, itp_vd, itp_def, itp_q, qav)
	""" Choose ap ≥ 0, bp ≥ 0 default with prob sr.gov[:repay] """
	jζ = 2 # IN REPAYMENT
	jdef = def_state(sr, jζ)
	jdef == false || throw(error("Wrong default state"))

	xguess = [guess[key][jζ] for key in [:b, :a]]
	xmin = [minimum(sr.gr[key]) for key in [:b, :a]]
	xmax = [maximum(sr.gr[key]) for key in [:b, :a]]
	xmax[1] = max(0.5*xmax[1], min(xmax[1], 2*state[:b]))

	# xguess *= 0.01
	c_guess, _ = budget_constraint_agg(sr, state, pz, pν, xguess, itp_def, itp_q, qav, jdef)
	if c_guess < 1e-2
		xguess[1] = state[:b] + 1e-2 * (.5*xmax[1]+.5*xmin[1] - state[:b])
		xguess[2] = 0.01
	end

	obj_f(x) = value(sr, state, pz, pν, x[1], x[2], itp_v, itp_vd, itp_def, itp_q, qav, jdef)
	# res = Optim.optimize(obj_f, xmin, xmax, xguess, Fminbox())
	# if !Optim.converged(res)
	# 	res = Optim.optimize(obj_f, xmin, xmax, xguess, Fminbox())
	# end
	# !Optim.converged(res) && println("WARNING: DIDN'T FIND SOL IN R")
	# x_opt = res.minimizer
	# bpv, apv = x_opt
	# vR = -obj_f(x_opt)

	# opt = Opt(:LN_COBYLA, length(xguess))
	# opt = Opt(:LD_MMA, length(xguess))
	opt = Opt(:LD_SLSQP, length(xguess))
	# opt = Opt(:LN_SBPLX, length(xguess))
	opt.lower_bounds = xmin
	opt.upper_bounds = xmax
	# opt.xtol_rel = 1e-16
	
	function F(x,g)
		if length(g) > 0
			g[:] = ForwardDiff.gradient(obj_f, x)
		end
		obj_f(x)
	end

	constr(x) = Euler_eq(sr, state, pz, pν, x[1], x[2], itp_ucT, itp_def, itp_q, qav, jdef)
	function G(x,g,v)
		if length(g) > 0
			g[:] = ForwardDiff.gradient(constr, x)
		end
		constr(x) - v
	end
	opt.max_objective = F
	opt.maxeval = 500
	if sr.opt[:Euler] == true
		inequality_constraint!(opt, (x,g) -> G(x,g,1e-6))
	end
	maxf, x_opt, ret = NLopt.optimize(opt, xguess)
	# println(ret)
	bpv, apv = x_opt
	vR = maxf

	ϕ = Dict(:a=>apv, :b=>bpv)

	ccv, cTv = budget_constraint_agg(sr, state, pz, pν, x_opt, itp_def, itp_q, qav, jdef)
	ϕ[:c] = ccv
	ϕ[:cT] = cTv
	return ϕ, vR
end

function opt_value_D(sr::SOEres, guess, state, pz, pν, itp_ucT, itp_v, itp_vd, itp_def, itp_q, qav)
	""" Choose ap between 0 and a, bp = bv, reenter mkts w prob θ """
	bpv = state[:b]
	jζ = 1 # IN DEFAULT
	jdef = def_state(sr, jζ)
	jdef == true || throw(error("Wrong default state"))

	xguess = [guess[key][jζ] for key in [:a]]
	xmin = [minimum(sr.gr[key]) for key in [:a]]
	xmax = [maximum(sr.gr[key]) for key in [:a]]

	c_guess, _ = budget_constraint_agg(sr, state, pz, pν, [bpv, xguess[1]], itp_def, itp_q, qav, jdef)
	if c_guess < 1e-2
		# xguess[1] = xmax[1] * 0.99
		xguess[1] = 1e-6
	end

	# obj_f(x) = -value(sr, state, pz, pν, bpv, x[1], itp_v, itp_vd, itp_def, itp_q, qav, jdef)
	# res = Optim.optimize(obj_f, xmin, xmax, xguess, Fminbox(GradientDescent()))
	obj_f(apv) = value(sr, state, pz, pν, bpv, first(apv), itp_v, itp_vd, itp_def, itp_q, qav, jdef)
	# res = Optim.optimize(obj_f, xmin, xmax, xguess, Fminbox())

	# res = Optim.optimize(obj_f, first(xmin), first(xmax), GoldenSection())

	# !Optim.converged(res) && println("WARNING: DIDN'T FIND SOL IN D")

	# opt = Opt(:LN_SBPLX, length(xguess))
	opt = Opt(:LD_SLSQP, length(xguess))
	opt.lower_bounds = xmin
	opt.upper_bounds = xmax
	
	function F(x,g)
		if length(g) > 0
			g[:] = ForwardDiff.gradient(obj_f, x)
		end
		obj_f(x)
	end

	constr(x) = Euler_eq(sr, state, pz, pν, bpv, first(x), itp_ucT, itp_def, itp_q, qav, jdef)
	function G(x,g,v)
		if length(g) > 0
			g[:] = ForwardDiff.gradient(constr, x)
		end
		constr(x) - v
	end
	opt.max_objective = F
	opt.maxeval = 500
	if sr.opt[:Euler] == true
		inequality_constraint!(opt, (x,g) -> G(x,g,1e-6))
	end
	maxf, x_opt, ret = NLopt.optimize(opt, xguess)
	# println(ret)
	apv = first(x_opt)
	vD = maxf

	ϕ = Dict(:a=>apv, :b=>bpv)
	
	ccv, cTv = budget_constraint_agg(sr, state, pz, pν, [bpv, apv], itp_def, itp_q, qav, jdef)
	ϕ[:c] = ccv
	ϕ[:cT] = cTv
	return ϕ, vD
end

function solve_optvalue!(new_v, new_ϕ, sr::SOEres, itp_ucT, itp_v, itp_vd, itp_def, itp_q)
	""" Loop over states and solve """
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
		ϕR, vR = opt_value_R(sr, guess, state, pz, pν, itp_ucT, itp_v, itp_vd, itp_def, itp_q, qaR)
		ϕD, vD = opt_value_D(sr, guess, state, pz, pν, itp_ucT, itp_v, itp_vd, itp_def, itp_q, qaD)

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
end

function prob_extreme_value(sr::SOEres, vR, vD)
	""" Apply extreme-value shocks formula """
	σV = sr.pars[:σV]
	prob = exp(vD/σV) / (exp(vR/σV) + exp(vD/σV))
	if isnan(prob)
		if vR > vD
			return 0.0
		else
			return 1.0
		end
	end
	return prob
end

function update_def!(sr::SOEres, new_v=sr.v)
	""" Computes default prob and value of entering period in repayment """
	""" Uses 'updated' sr.v[:D] and sr.v[:R] by default """
	itp_vd = make_itp(sr,new_v[:D]);
	ℏ = sr.pars[:ℏ]

	Jgrid = agg_grid(sr);
	for js in 1:size(Jgrid,1)
		jv = Jgrid[js,:]

		state = S(sr, jv)
		st = [state[key] for key in statenames(sr)]
		st_def = corr(sr, st)

		vR = new_v[:R][jv...]
		vD = itp_vd(st_def...)

		def_prob = prob_extreme_value(sr,vR,vD)

		new_v[:def][jv...] = def_prob
		new_v[:V][jv...] = def_prob * vD + (1-def_prob) * vR
	end
end

function vfi_iter!(new_v, new_ϕ, sr::SOEres)
	""" Interpolate values and prices to use as next period values """
	itp_v   = make_itp(sr, sr.v[:V]);
	itp_vd  = make_itp(sr, sr.v[:D]);
	itp_def = make_itp(sr, sr.v[:def]);
	itp_q   = make_itp(sr, sr.eq[:qb]);
	itp_ucT = make_itp(sr, sr.ϕ[:cT]);

	solve_optvalue!(new_v, new_ϕ, sr, itp_ucT, itp_v, itp_vd, itp_def, itp_q);
	update_def!(sr, new_v)
	
	nothing
end

function update_sr!(y, new_y, upd_η = 1)
	for key in keys(y)
		y[key] = y[key] + upd_η * (new_y[key] - y[key])
	end
	nothing
end

function vfi!(sr::SOEres; tol::Float64=5e-4, maxiter::Int64=1000, verbose::Bool=false)
	""" Main Loop """
	iter, dist = 0, 1+tol
	avg_time = 0.0
	dist_v, dist_ϕ = zeros(2)

	new_v = Dict(key => similar(val) for (key, val) in sr.v)
	new_ϕ = Dict(key => similar(val) for (key, val) in sr.ϕ)

	upd_η = 0.9
	t0 = time()
	while dist > tol && iter < maxiter
		iter += 1
		verbose && print("Iteration $iter\n")

		old_q = copy(sr.eq[:qb]);
		""" Update debt prices (for use as next period prices) """
		update_q!(sr, verbose = verbose)
		dist_q = norm(sr.eq[:qb]-old_q) / max(1,norm(old_q))
		# sr.eq[:qb] = old_q + upd_ηq * (sr.eq[:qb] - old_q)
		norm_q = norm(sr.eq[:qb])

		t1 = time()
		""" Iterate on the value functions """
		vfi_iter!(new_v, new_ϕ, sr);
		t = time() - t1

		avg_time = (avg_time * (iter - 1) + t) / iter

		dist_v = maximum([ norm(new_v[key] - sr.v[key]) / max(1,norm(sr.v[key])) for key in keys(sr.v) ])
		dist_ϕ = maximum([ norm(new_ϕ[key] - sr.ϕ[key]) / max(1,norm(sr.ϕ[key])) for key in keys(sr.ϕ) ])

		norm_v = maximum([norm(sr.v[key]) for key in keys(sr.v)])
		norm_ϕ = maximum([norm(sr.ϕ[key]) for key in keys(sr.ϕ)])
		
		dist = max(dist_v, dist_ϕ, dist_q)

		for key in keys(sr.ϕ)
			print("$key: $(norm(new_ϕ[key] - sr.ϕ[key]) / max(1,norm(sr.ϕ[key]))) \n")
		end
		# for (key, val) in sr.ϕ
		# 	print("||$(key)|| = $(norm(val))\n")
		# end

		update_sr!(sr.v, new_v, upd_η)
		update_def!(sr) # To update [:def] and [:V] consistently with [:D], [:R]
		update_sr!(sr.ϕ, new_ϕ, upd_η)

		upd_η = max(upd_η * 0.995, 0.4)


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
