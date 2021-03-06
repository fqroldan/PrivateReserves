using Printf, Interpolations

function repayment_def(sr::SOEres, qv, qvd, jζ, jζp)
	""" Use future q and reentry/default probs to calculate repayment """
	θ, κ, δ, ℏ = [sr.pars[sym] for sym in [:θ, :κ, :δ, :ℏ]]

	rep = 0.0
	if def_state(sr, jζp) == false # Reentry to markets or normal repayment
		rep = κ + (1-δ) * qv
	elseif def_state(sr, jζ) == true # Default in t and t+1
		rep = qv
	elseif def_state(sr, jζ) == false # Defaulted at t+1
		rep = (1-ℏ) * qvd
	end
	return rep
end

function prob_def(sr::SOEres, jζ, jζp, defp)
	""" Use repayment policy and reentry params to calculate prob of default states """
	def_prob = 1.0
	if def_state(sr, jζ) == false # Starting from normal access (ζv == 1)
		if def_state(sr, jζp) == true # Default probability
			def_prob = defp
		else # Repayment probability
			def_prob = 1-defp
		end
	else # Starting from default
		if def_state(sr, jζp) == true # Remain probability
			def_prob = 1-sr.pars[:θ]
		else # Reentry probability
			def_prob = sr.pars[:θ]
		end
	end
	return def_prob
end

function iterate_q!(new_qa, new_qb, sr::SOEres, itp_q, itp_def)
	""" One iteration of the debt price """
	ψ, σz, r, ℏ = [sr.pars[key] for key in [:ψ, :σz, :r, :ℏ]]

	Jgrid = agg_grid(sr)
	Threads.@threads for js in 1:size(Jgrid,1)
		jv = Jgrid[js, :]
		jd = S_index(sr, jv)
		
		state = S(sr, jv)
		zv = state[:z]
		νv = state[:ν]

		pz = sr.prob[:z][jd[:z],:]
		pν = sr.prob[:ν][jd[:ν],:]

		for jζ in 1:2
			bpv = sr.ϕ[:b][jv...,jζ]
			apv = sr.ϕ[:a][jv...,jζ]

			qna = 0.0
			qnb = 0.0
			for (jzp, zpv) in enumerate(sr.gr[:z]), (jνp, νpv) in enumerate(sr.gr[:ν])
				defp = itp_def(bpv, apv, zpv, νpv)
				prob = pz[jzp] * pν[jνp]

				ϵpv = innov_z(sr, zpv, zv)
				sdf = SDF(sr, νv, ϵpv)

				for jζp in 1:2
					ζpv = sr.gr[:def][jζp]
					qvp = itp_q(bpv, apv, zpv, νpv, ζpv)
					qvd = itp_q((1-ℏ)*bpv, apv, zpv, νpv, ζpv)
					qnb += prob * prob_def(sr, jζ, jζp, defp) * repayment_def(sr, qvp, qvd, jζ, jζp) * sdf
				end
				qna += prob * sdf * 1
			end
			new_qa[jv..., jζ] = qna
			new_qb[jv..., jζ] = qnb
		end
	end
end

function update_q!(sr::SOEres; tol::Float64=1e-6, maxiter::Int64=500, verbose::Bool=true)
	""" Iterates on the debt price """
	iter, dist = 0, 1+tol
	upd_η = 0.75

	new_qa = zeros(size(sr.eq[:qb]))
	new_qb = zeros(size(sr.eq[:qb]))

	t0 = time()
	while dist > tol && iter < maxiter
		iter += 1

		old_q   = copy(sr.eq[:qb])
		itp_q   = make_itp(sr, sr.eq[:qb]);
		itp_def = make_itp(sr, sr.v[:def]);

		iterate_q!(new_qa, new_qb, sr, itp_q, itp_def)

		dist = sqrt( sum((new_qb - old_q).^2) / sum(old_q.^2) )

		sr.eq[:qb] = old_q + upd_η * (new_qb - old_q)
		sr.eq[:qa] = new_qa
	end

	# sr.eq[:qa] = ones(size(sr.eq[:qa])) * exp(-sr.pars[:r])

	tT = time()
	if dist <= tol
		if verbose
			print("Updated prices after $iter iterations in $(time_print(tT-t0))\n")
		end
	else
		print("WARNING: Iteration on qᵍ aborted at distance $(@sprintf("%.3g",dist)) after $(time_print(tT-t0))\n")
	end

	nothing
end

function get_eqm(sr::SOEres, bp, ap, state, pz, pν, jζ, itp_def, qav, itp_q)
	ζv = sr.gr[:def][jζ]
	jdef = def_state(sr, jζ)

	cT = budget_constraint_T(sr, state, pz, pν, [bp, ap], itp_def, itp_q, qav, jdef)
	cT = max(0.0, cT)
	hp = eq_h(sr, cT, jdef)
	yN = prod_N(sr, hp, jdef)
	yT = output_T(sr, state, jdef)

	output = CES_aggregator(sr, yT, yN)
	CA = yT - cT

	new_eq = Dict([(:cT,cT), (:cN,yN), (:output,output), (:labor,hp), (:CA,CA)])
end

function update_eqm!(sr::SOEres)
	itp_def = make_itp(sr, sr.v[:def])
	itp_q = make_itp(sr, sr.eq[:qb])

	Jgrid = agg_grid(sr)
	Threads.@threads for js in 1:size(Jgrid,1)
		jv = Jgrid[js, :]
		jd = S_index(sr, jv)
		
		state = S(sr, jv)
		zv = state[:z]
		νv = state[:ν]
		pz = sr.prob[:z][jd[:z],:]
		pν = sr.prob[:ν][jd[:ν],:]

		for jζ in 1:2
			bp, ap = [sr.ϕ[key][jv..., jζ] for key in [:b, :a]]

			qav = 0.0
			if def_state(sr, jζ) == 1
				qav = sr.eq[:qa][jv..., jζ]
			else
				qav = sr.eq[:qa][jv..., jζ]
			end

			new_eq = get_eqm(sr, bp, ap, state, pz, pν, jζ, itp_def, qav, itp_q)
			for (key, val) in new_eq
				if haskey(sr.eq, key)
					sr.eq[key][jv..., jζ] = val
				end
			end
		end
	end
end




