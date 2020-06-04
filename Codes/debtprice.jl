using Printf, Interpolations

function repayment_def(sr::SOEres, qv, jζ, jζp)
	""" Use future q and reentry/default probs to calculate repayment """
	θ, κ, δ, ℏ = [sr.pars[sym] for sym in [:θ, :κ, :δ, :ℏ]]

	rep = 0.0
	if jζp == 2 # Reentry to markets or normal repayment
		rep = κ + (1-δ) * qv
	elseif jζ == 1 # Default in t and t+1
		rep = qv
	elseif jζ == 2 # Defaulted at t+1
		rep = (1-ℏ) * qv
	end
	return rep
end

function prob_def(sr::SOEres, jζ, jζp, defp)
	""" Use repayment policy and reentry params to calculate prob of default states """
	def_prob = 1.0
	if jζ == 2 # Starting from normal access (ζv == 1)
		def_state(sr, jζ) == false || throw(error("Wrong def"))
		if jζp == 1 # Default probability
			def_state(sr, jζp) == true || throw(error("Wrong def"))
			def_prob = defp
		else # Repayment probability
			def_state(sr, jζp) == false || throw(error("Wrong def"))
			def_prob = 1-defp
		end
	else # Starting from default
		if jζp == 1 # Remain probability
			def_state(sr, jζp) == true || throw(error("Wrong def"))
			def_prob = 1-sr.pars[:θ]
		else # Reentry probability
			def_state(sr, jζp) == false || throw(error("Wrong def"))
			def_prob = sr.pars[:θ]
		end
	end
	return def_prob
end

function iterate_q(sr::SOEres, itp_q, itp_def)
	""" One iteration of the debt price """
	new_q = zeros(size(sr.eq[:qb]))
	ψ, σz, r = [sr.pars[key] for key in [:ψ, :σz, :r]]

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

			qn = 0.0
			for (jzp, zpv) in enumerate(sr.gr[:z]), (jνp, νpv) in enumerate(sr.gr[:ν])
				defp = itp_def(bpv, apv, zpv, νpv)
				prob = pz[jzp] * pν[jνp]

				ϵpv = innov_z(sr, zpv, zv)
				sdf = exp(-r - νv * (ψ * ϵpv + 0.5 * ψ^2*σz^2))

				for jζp in 1:2
					ζpv = sr.gr[:def][jζp]
					qvp = itp_q(bpv, apv, zpv, νpv, ζpv)
					qn += prob * prob_def(sr, jζ, jζp, defp) * repayment_def(sr,qvp, jζ, jζp) * sdf
				end
			end
			new_q[jv..., jζ] = qn
		end
	end
	return new_q
end

function update_q!(sr::SOEres; tol::Float64=1e-6, maxiter::Int64=500, verbose::Bool=true)
	""" Iterates on the debt price """
	iter, dist = 0, 1+tol
	upd_η = 0.75

	t0 = time()
	while dist > tol && iter < maxiter
		iter += 1

		old_q   = copy(sr.eq[:qb])
		itp_q   = make_itp(sr, sr.eq[:qb]);
		itp_def = make_itp(sr, sr.v[:def]);

		new_q = iterate_q(sr, itp_q, itp_def)

		dist = sqrt( sum((new_q - old_q).^2) / sum(old_q.^2) )

		sr.eq[:qb] = old_q + upd_η * (new_q - old_q)
	end

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


function update_rep!(sr::SOEres)
	""" Updates the repayment probability given (S_t, z_{t+1}, ν_{t+1}) """
	itp_def = make_itp(sr, sr.v[:def])

	Jgrid = agg_grid(sr)
	jζ = 2 # Default choice at t+1 only if normal access at t
	for jj in 1:size(Jgrid,1)
		js = Jgrid[jj,:]
		bpv, apv = [sr.ϕ[key][js...,jζ] for key in [:b, :a]]
		for (jzp, zpv) in enumerate(sr.gr[:z]), (jνp, νpv) in enumerate(sr.gr[:ν])
			def_prob = max(0,min(1,itp_def(bpv, apv, zpv, νpv)))
			sr.gov[:repay][js...,jzp,jνp] = 1 - def_prob
		end
	end
	nothing
end
