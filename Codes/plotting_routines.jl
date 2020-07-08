""" Define styles """

def_style = let
	axis = attr(showgrid = true, gridcolor="#e2e2e2", gridwidth=0.5, zeroline=false)
	layout = Layout(xaxis = axis, yaxis=axis)
	Style(layout=layout)
end

slides_def = let
	layout = Layout(plot_bgcolor="#fafafa", paper_bgcolor="#fafafa",
		width=1920*0.45, height=1080*0.45, font_size=16, font_family="Lato",
		legend = attr(orientation = "h", x=0.05))
	Style(def_style, layout=layout)
end

dark_bg = let
	axis = attr(gridcolor="#1b1b1b")
	layout = Layout(plot_bgcolor="#020202", paper_bgcolor="#020202", font_color="white", xaxis=axis,yaxis=axis)
	Style(layout=layout)
end
slides_dark = Style(slides_def, dark_bg)

paper = let
	layout = Layout(width = 1920 * 0.5, height = 1080 * 0.35, font_size=16, font_family = "Linux Libertine",
		legend = attr(orientation = "h", x=0.05))
	Style(def_style, layout=layout)
end

contsty(;diverg::Bool=false) = let
	colpal = ifelse(diverg, ColorSchemes.davos, ColorSchemes.lajolla)
	colscale = [[vv, get(colpal, vv)] for vv in range(0,1,length=100)]
	c_att = attr(colorscale=colscale, autocontour=false)
	Style(trace=Dict(:contour=>c_att))
end

default_eval_vec(sr::SOEres) = floor(Int, N(sr,:b)*0.5), floor(Int, N(sr,:a)*0.5), floor(Int, N(sr,:z)*0.5), floor(Int, N(sr,:ν)*0.5)

function default_eval_points(sr::SOEres)
	jb, ja, jz, jν = default_eval_vec(sr)
	return Dict(:b=>jb, :a=>ja, :z=>jz, :ν=>jν)
end

function grab_vec(sr::SOEres, y::Array{Float64,K}, key::Symbol, jζ=2; eval_points::Dict{Symbol, Int64}=Dict{Symbol,Int64}()) where K
	jk = findfirst(statenames(sr).==key)

	eval_indices = default_eval_points(sr)
	for (key, val) in eval_points
		if haskey(eval_indices, key)
			eval_indices[key] = val
		end
	end

	jvdef = [eval_indices[key] for key in statenames(sr)]
	if length(jvdef) + 1 == K
		jvdef = push!(jvdef, jζ)
	end

	yv = zeros(size(y,jk))
	for jj in 1:size(y,jk)
		jv = jvdef
		jv[jk] = jj

		yv[jj] = y[jv...]
	end
	return yv
end

function grab_mat(sr::SOEres, y::Array{Float64,K}, k1::Symbol, k2::Symbol, jζ=2; eval_points::Dict{Symbol, Int64}=Dict{Symbol,Int64}()) where K
	j1 = findfirst(statenames(sr).==k1)
	N1 = size(y,j1)
	j2 = findfirst(statenames(sr).==k2)
	N2 = size(y,j2)

	eval_indices = default_eval_points(sr)
	for (key, val) in eval_points
		if haskey(eval_indices, key)
			eval_indices[key] = val
		end
	end

	jvdef = [eval_indices[key] for key in statenames(sr)]
	if length(jvdef) + 1 == K
		jvdef = push!(jvdef, jζ)
	end

	ym = zeros(N1, N2)
	for jx in 1:N1, jy in 1:N2
		jv = jvdef
		jv[j1] = jx
		jv[j2] = jy

		ym[jx, jy] = y[jv...]
	end
	return ym
end

makeplot_eq(sr::SOEres, ykey::Vector{Symbol}, xkey::Symbol, def::Bool=true; style::Style=slides_def) = makeplot_sr(sr, sr.eq, ykey, xkey, def; style=style)

makeplot_v(sr::SOEres, ykey::Vector{Symbol}, xkey::Symbol, def::Bool=false; style::Style=slides_def) = makeplot_sr(sr, sr.v, ykey, xkey, def; style=style)

makeplot_sr(sr::SOEres, srdict::Dict{Symbol, Array{Float64, K}}, ykey::Symbol, xkey::Symbol, def::Bool=true; style::Style=slides_def, eval_points::Dict{Symbol,Int64}=Dict{Symbol,Int64}()) where K = makeplot_sr(sr, srdict, [ykey], xkey, def; style=slides_def, eval_points=eval_points)

makeplot_eq(sr::SOEres, ykey::Symbol, xkey::Symbol, def::Bool=true; style::Style=slides_def) = makeplot_eq(sr, [ykey], xkey, def; style=style)

makeplot_v(sr::SOEres, ykey::Symbol, xkey::Symbol, def::Bool=true; style::Style=slides_def) = makeplot_v(sr, [ykey], xkey, def; style=style)

function makeplot_sr(sr::SOEres, srdict::Dict{Symbol, Array{Float64, K}}, ykey::Vector{Symbol}, xkey::Symbol; def::Bool=true, style::Style=slides_def, eval_points::Dict{Symbol,Int64}=Dict{Symbol,Int64}()) where K

	Nζ = 1+def

	xvec = sr.gr[xkey]
	yvec = [[grab_vec(sr, srdict[y], xkey, jj, eval_points=eval_points) for jj in 1:Nζ] for (jy,y) in enumerate(ykey)]

	namevec = [[string(y) * ifelse(Nζ==2,ifelse(def_state(sr, jζ), ": Default", ": Repayment"),"") for jζ in 1:2] for (jy,y) in enumerate(ykey)]

	data = [scatter(x=xvec, y=yvec[jy][jζ], name = namevec[jy][jζ]) for jζ in Nζ:-1:1, jy in 1:length(ykey)]

	p1 = plot(data[:],style=style, Layout(title=uppercasefirst(string(first(ykey))), xaxis_title="<i>"*string(xkey)))

	return p1
end

function make_contour(sr::SOEres, srd::Dict, zkey::Symbol, xkey::Symbol, ykey::Symbol; style::Style=slides_def, eval_points::Dict{Symbol,Int64}=Dict{Symbol,Int64}())
	y = grab_mat(sr, srd[zkey], xkey, ykey, eval_points=eval_points)

	data = contour(x=sr.gr[xkey], y=sr.gr[ykey], z=y)

	layout = Layout(xaxis_title="<i>"*string(xkey), yaxis_title="<i>"*string(ykey))

	return plot(data, style=Style(style, contsty(diverg=false)), layout)
end

function make_comp_V(sr::SOEres, xkey::Symbol, ykey::Symbol=:nothing; style::Style=slides_def, eval_points::Dict{Symbol,Int64}=Dict{Symbol,Int64}())
	ℏ = sr.pars[:ℏ]
	index_b = findfirst(statenames(sr).==:b)

	itp_vd = make_itp(sr,sr.v[:D]);
	new_v = Dict(key => copy(val) for (key, val) in sr.v)

	Jgrid = agg_grid(sr);
	for js in 1:size(Jgrid,1)
		jv = Jgrid[js,:]

		state = [S(sr, jv)[key] for key in statenames(sr)]
		corr = ones(length(state))
		corr[index_b] *= (1-ℏ)
		eval_state = state .* corr

		vR = new_v[:R][jv...]
		vD = itp_vd(eval_state...)

		def_prob = prob_extreme_value(sr,vR,vD)

		new_v[:def][jv...] = def_prob
		new_v[:V][jv...] = def_prob * vD + (1-def_prob) * vR
		new_v[:R][jv...] = vR
		new_v[:D][jv...] = vD
	end
	
	if ykey == :nothing
		makeplot_sr(sr, new_v, [:R, :D], xkey, def=false, style=style, eval_points=eval_points)
	else
		yR = grab_mat(sr, new_v[:R], xkey, ykey, eval_points=eval_points)
		yD = grab_mat(sr, new_v[:D], xkey, ykey, eval_points=eval_points)

		maxz = maximum(abs,extrema(yR-yD))

		data = contour(x=sr.gr[xkey], y=sr.gr[ykey], z=yR-yD,
			# contours=Dict(:start=>-maxz, :end=>maxz)
			)
		layout = Layout(xaxis_title="<i>"*string(xkey), yaxis_title="<i>"*string(ykey), title="<i>V<sup>R</sup> - V<sup>D</sup>")
		return plot(data, style=Style(style, contsty(diverg=true)), layout)
	end
end

make_debtprice(sr::SOEres, xkey::Symbol, ykey::Symbol; style::Style=slides_def, eval_points::Dict{Symbol,Int64}=Dict{Symbol,Int64}()) = make_contour(sr, sr.eq, :qb, xkey, ykey, style=style, eval_points=eval_points)

make_defprob(sr::SOEres, xkey::Symbol, ykey::Symbol; style::Style=slides_def, eval_points::Dict{Symbol,Int64}=Dict{Symbol,Int64}()) = make_contour(sr, sr.v, :def, xkey, ykey, style=style, eval_points=eval_points)


function plot_simul(sr::SOEres, pp::Path)
	defs = [tt for tt in 1:periods(pp) if series(pp, :newdef)[tt] == 1 || (tt == 1 && series(pp, :ζ)[1] == 0)];
	recs = [tt for tt in 1:periods(pp) if series(pp, :reentry)[tt] == 1 || (tt == periods(pp) && series(pp,:ζ)[end] == 0)];

	print("$(length(defs)) defaults in $(floor(Int,periods(pp)/4)) years.\n")

	shapes = [attr(type="rect", xref = "x", x0 = defs[jj]./4, x1=recs[jj]./4, yref ="paper", y0=0, y1=1, fillcolor="#bababa", opacity=0.2, line_width=0) for jj in 1:length(defs)];

	data = [
		scatter(x=(1:periods(pp))./4, yaxis="y3", y=100series(pp,:a)./(4*series(pp,:output)), name="<i>a/y")
		scatter(x=(1:periods(pp))./4, yaxis="y3", y=100series(pp,:b)./(4*series(pp,:output)), name="<i>b/y")
		scatter(x=(1:periods(pp))./4, yaxis="y2", y=100*(1 .- series(pp,:labor)), name="<i>unemp")
		scatter(x=(1:periods(pp))./4, yaxis="y1", y=exp.(series(pp,:z)), name="<i>TFP")
		]
	layout = Layout(shapes=shapes,
		yaxis1 = attr(domain=[0, 0.3]),
		yaxis2 = attr(domain=[0.33, 0.67], title="<i>%"),
		yaxis3 = attr(domain=[0.7, 1], title="<i>% of GDP"),
		)

	plot(data, layout, style=slides_def)
end