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

function price_debt(sr::SOEres, xp, zv, νv, pz, pν, itp_def, itp_q)
	""" Iterates once on the debt price using next period's state """
	δ, ℏ, ψ, σz, r = [sr.pars[sym] for sym in [:δ, :ℏ, :ψ, :σz, :r]]
	qb = 0.0
	bp, ap = xp
	for (jzp, zpv) in enumerate(sr.gr[:z]), (jνp, νpv) in enumerate(sr.gr[:ν])
		prob = pz[jzp] * pν[jνp]

		ϵpv = innov_z(sr, zpv, zv)
		sdf = exp(-r - νv * (ψ * ϵpv - 0.5 * ψ^2*σz^2))

		jζp = 1 # Default
		rep_default= (1-δ) * (1-ℏ) * itp_q(bp*(1-ℏ),ap,zpv,νpv, sr.gr[:def][jζp])
		jζp = 2 # Repayment
		rep_normal = δ + (1-δ) * itp_q(bp,ap,zpv,νpv, sr.gr[:def][jζp])

		def_prob = itp_def(bp,ap,zpv,νpv)

		qb += prob * sdf * (def_prob * rep_default + (1-def_prob) * rep_normal)
	end
	return qb
end

function budget_constraint_T(sr::SOEres, state, pz, pν, xp, itp_def, itp_q, jdef::Bool)
	""" Computes consumption of T given state and choices of new debt and reserves """

	yT = output_T(sr, state, jdef)
	δv = bond_decay(sr, jdef)

	bv, av, zv, νv = [state[key] for key in [:b, :a, :z, :ν]]
	qa = exp(-sr.pars[:r])

	bp, ap = xp

	debt_operations = 0.0
	if !jdef
		qb = price_debt(sr, xp, zv, νv, pz, pν, itp_def, itp_q)
		debt_operations = qb * (bp - (1-δv) * bv) - δv * bv
	end

	cT = yT + av - qa*ap + debt_operations
	return cT
end

function budget_constraint_agg(sr::SOEres, state, pz, pν, xp, itp_def, itp_q, jdef::Bool)
	""" Computes aggregate of C at state and choices of new debt and reserves """
	cT = budget_constraint_T(sr, state, pz, pν, xp, itp_def, itp_q, jdef)

	cT = max(0.0, cT)

	h = eq_h(sr, cT)
	yN = prod_N(sr, h)
	c = CES_aggregator(sr, cT, yN)

	return c
end

function value(sr::SOEres, state, pz, pν, bp, ap, itp_v, itp_vd, itp_def, itp_q, jdef::Bool)
	""" Computes V given state and choices of consumption, debt, reserves """
	θ, ℏ, β = [sr.pars[sym] for sym in [:θ, :ℏ, :β]]

	# xp = Dict([(:b,bp), (:a,ap)])

	xp = [bp, ap]

	# bp = max(bp, minimum(sr.gr[:b]))
	# bp = min(bp, maximum(sr.gr[:b]))

	c = budget_constraint_agg(sr, state, pz, pν, xp, itp_def, itp_q, jdef)
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

function opt_value_R(sr::SOEres, guess, state, pz, pν, itp_v, itp_vd, itp_def, itp_q)
	""" Choose ap ≥ 0, bp ≥ 0 default with prob sr.gov[:repay] """
	jζ = 2 # IN REPAYMENT
	jdef = def_state(sr, jζ)
	jdef == false || throw(error("Wrong default state"))

	xguess = [guess[key][jζ] for key in [:b, :a]]
	xmin = [minimum(sr.gr[key]) for key in [:b, :a]]
	xmax = [maximum(sr.gr[key]) for key in [:b, :a]]

	obj_f(x) = -value(sr, state, pz, pν, x[1], x[2], itp_v, itp_vd, itp_def, itp_q, jdef)
	res = Optim.optimize(obj_f, xmin, xmax, xguess, Fminbox(NelderMead()))

	if !Optim.converged(res)
		res = Optim.optimize(obj_f, xmin, xmax, xguess, Fminbox(BFGS()))
	end
	
	!Optim.converged(res) && println("WARNING: DIDN'T FIND SOL IN R")
	
	x_opt = res.minimizer
	bpv, apv = x_opt
	vR = -obj_f(x_opt)

	ϕ = Dict(:a=>apv, :b=>bpv)

	ccv = budget_constraint_agg(sr, state, pz, pν, x_opt, itp_def, itp_q, jdef)
	ϕ[:c] = ccv
	return ϕ, vR
end

function opt_value_D(sr::SOEres, guess, state, pz, pν, itp_v, itp_vd, itp_def, itp_q)
	""" Choose ap between 0 and a, bp = bv, reenter mkts w prob θ """
	bpv = state[:b]
	jζ = 1 # IN DEFAULT
	jdef = def_state(sr, jζ)
	jdef == true || throw(error("Wrong default state"))

	xguess = [guess[key][jζ] for key in [:a]]
	xmin = [minimum(sr.gr[key]) for key in [:a]]
	xmax = [maximum(sr.gr[key]) for key in [:a]]

	obj_f(x) = -value(sr, state, pz, pν, bpv, x[1], itp_v, itp_vd, itp_def, itp_q, jdef)
	res = Optim.optimize(obj_f, xmin, xmax, xguess, Fminbox(GradientDescent()))

	!Optim.converged(res) && println("WARNING: DIDN'T FIND SOL IN D")

	x_opt = res.minimizer
	apv = first(x_opt)
	vD = -obj_f(x_opt)

	ϕ = Dict(:a=>apv, :b=>bpv)
	
	ccv = budget_constraint_agg(sr, state, pz, pν, [bpv, apv], itp_def, itp_q, jdef)
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

		""" New xp and values in repayment and default """
		ϕR, vR = opt_value_R(sr, guess, state, pz, pν, itp_v, itp_vd, itp_def, itp_q)
		ϕD, vD = opt_value_D(sr, guess, state, pz, pν, itp_v, itp_vd, itp_def, itp_q)

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
	κ = sr.pars[:κ]
	return exp(vD/κ) / (exp(vR/κ) + exp(vD/κ))
end

function update_def!(sr::SOEres, new_v)
	""" Computes default prob and value of entering period in repayment """
	itp_vd = make_itp(sr,new_v[:D]);
	ℏ = sr.pars[:ℏ]

	Jgrid = agg_grid(sr);
	for js in 1:size(Jgrid,1)
		jv = Jgrid[js,:]

		state = S(sr, jv)
		st = [state[key] for key in [:b,:a,:z,:ν]]
		st_def = corr(sr, st)
		# corr = ones(length(st))
		# index_b = findfirst(statenames(sr).==:b)
		# corr[index_b] *= (1-ℏ)
		# st_def = st .* corr

		vR = new_v[:R][jv...]
		vD = itp_vd(st_def...)

		def_prob = prob_extreme_value(sr,vR,vD)

		new_v[:def][jv...] = def_prob
		new_v[:V][jv...] = def_prob * vD + (1-def_prob) * vR
	end
end

function vfi_iter(sr::SOEres)
	""" Interpolate values and prices to use as next period values """
	itp_v  = make_itp(sr, sr.v[:V]);
	itp_vd = make_itp(sr, sr.v[:D]);
	itp_def = make_itp(sr, sr.v[:def]);
	itp_q  = make_itp(sr, sr.eq[:qb]);

	new_v, new_ϕ = solve_optvalue(sr, itp_v, itp_vd, itp_def, itp_q);
	update_def!(sr, new_v)

	return new_v, new_ϕ
end

function update_sr!(y, new_y)
	upd_η = 0.75
	for key in keys(y)
		y[key] = y[key] + upd_η * (new_y[key] - y[key])
	end
	nothing
end

function vfi!(sr::SOEres; tol::Float64=1e-4, maxiter::Int64=500, verbose::Bool=false)
	""" Main Loop """
	iter, dist = 0, 1+tol
	avg_time = 0.0
	dist_v, dist_ϕ = zeros(2)

	t0 = time()
	while dist > tol && iter < maxiter
		iter += 1

		old_q = copy(sr.eq[:qb])
		""" Update debt prices (for use as next period prices) """
		update_q!(sr, verbose = false)
		dist_q = sum( (sr.eq[:qb]-old_q).^2 ) / sum(old_q.^2)

		old_v = copy(sr.v)
		old_ϕ = copy(sr.ϕ)

		t1 = time()
		""" Iterate on the value functions """
		new_v, new_ϕ = vfi_iter(sr)
		t = time() - t1

		avg_time = (avg_time * (iter - 1) + t) / iter

		dist_v = maximum([ sum((new_v[key] - old_v[key]).^2) / sum(old_v[key].^2) for key in keys(sr.v) ])
		dist_ϕ = maximum([ sum((new_ϕ[key] - old_ϕ[key]).^2) / sum(old_ϕ[key].^2) for key in keys(sr.ϕ) ])

		dist = max(dist_v, dist_ϕ, dist_q)

		update_sr!(sr.v, new_v)
		update_sr!(sr.ϕ, new_ϕ)

		# for key in keys(sr.v)
		# 	print("||$(key)|| = $(norm(sr.v[key]))\n")
		# end

		if verbose && iter % 10 == 0
			print("After $iter iterations (avg time = $(time_print(avg_time))), d(v,ϕ,q) = $(@sprintf("%0.3g",dist_v)), $(@sprintf("%0.3g",dist_ϕ)), $(@sprintf("%0.3g",dist_q)) \n")
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
