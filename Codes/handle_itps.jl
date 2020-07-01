using Interpolations

const itp_sr{K} = Interpolations.GriddedInterpolation{Float64,K,Float64,Gridded{Linear},NTuple{K,Array{Float64,1}}}

function make_itp(sr::SOEres{K,K1}, y::Array) where {K,K1}

	if length(size(y)) == K
		knots = (sr.gr[:b], sr.gr[:a], sr.gr[:z], sr.gr[:ν])
	elseif length(size(y)) == K1
		knots = (sr.gr[:b], sr.gr[:a], sr.gr[:z], sr.gr[:ν], sr.gr[:def])
	# elseif length(size(y)) == K2
	# 	knots = (sr.gr[:b], sr.gr[:a], sr.gr[:z], sr.gr[:ν], sr.gr[:z], sr.gr[:ν])
	end

	itp = interpolate(knots, y, Gridded(Linear()))

	itp_obj = extrapolate(itp, Interpolations.Line())

	return itp_obj
end