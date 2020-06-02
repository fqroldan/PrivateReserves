using Printf, Interpolations

function repayment_def(sr::SOEres, qv, jζ, jζp)
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

function prob_def(sr::SOEres, jζ, jζp, repv)
	rep_prob = 1.0
	if jζ == 2 # Starting from normal access (ζv == 1)
		if jζp == 1 # Default probability
			rep_prob = 1-repv
		else # Repayment probability
			rep_prob = repv
		end
	else # Starting from default
		if jζp == 1 # Remain probability
			rep_prob = 1-sr.pars[:θ]
		else # Reentry probability
			rep_prob = sr.pars[:θ]
		end
	end
	return rep_prob
end

function iterate_q(sr::SOEres, itp_q::itp_sr)
	new_q = zeros(size(sr.eq[:qb]))

	Jgrid = agg_grid(sr)
	Threads.@threads for js in 1:size(Jgrid,1)
		jb, ja, jz, jϵ = Jgrid[js, :]

		pz = sr.prob[:z][jz,:]
		pϵ = sr.prob[:ϵ][jϵ,:]

		bpv_vec = [sr.ϕ[:b][jb,ja,jz,jϵ,jj] for jj in 1:2]
		apv_vec = [sr.ϕ[:a][jb,ja,jz,jϵ,jj] for jj in 1:2]
		rep = zeros(2)
		for (jzp, zpv) in enumerate(sr.gr[:z]), (jϵp, ϵpv) in enumerate(sr.gr[:ϵ])
			prob = pz[jzp] * pϵ[jϵp]
			repv = sr.gov[:repay][jb,ja,jz,jϵ,jzp,jϵp]

			for jζp in 1:2, jζ in 1:2
				ζpv = sr.gr[:def][jζp]
				qv = itp_q(bpv_vec[jζp], apv_vec[jζp], zpv, ϵpv, ζpv)

				rep[jζ] += prob * prob_def(sr, jζ, jζp, repv) * repayment_def(sr, qv, jζ, jζp)
			end
		end

		new_q[jb, ja, jz, jϵ, :] .= rep
	end
	return new_q
end

function update_q!(sr::SOEres; tol::Float64=1e-6, maxiter::Int64=500, verbose::Bool=true)
	iter, dist = 0, 1+tol
	upd_η = 0.75

	t0 = time()
	while dist > tol && iter < maxiter
		iter += 1

		old_q = copy(sr.eq[:qb])
		itp_q = make_itp(sr, sr.eq[:qb]);

		new_q = iterate_q(sr, itp_q)

		dist = sqrt( sum((new_q - old_q).^2) / sum(old_q.^2) )

		sr.eq[:qb] = old_q + upd_η * (new_q - old_q)
	end

	tT = time()
	if dist <= tol
		if verbose
			print_save("Updated prices after $iter iterations in $(time_print(tT-t0))")
		end
	else
		print_save("WARNING: Iteration on qᵍ aborted at distance $(@sprintf("%.3g",dist)) after $(time_print(tT-t0))")
	end

	nothing
end
