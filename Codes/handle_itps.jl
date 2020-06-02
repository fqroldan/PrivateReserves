using Interpolations

const itp_sr{K} =  Interpolations.GriddedInterpolation{Float64,K,Float64,Gridded{Linear},NTuple{K,Array{Float64,1}}}

function make_itp(sr::SOEres{K,K1,K2}, y::Array) where {K,K1,K2}

	if length(size(y)) == K
		knots = (sr.gr[:b], sr.gr[:a], sr.gr[:z], sr.gr[:ϵ])
	elseif length(size(y)) == K1
		knots = (sr.gr[:b], sr.gr[:a], sr.gr[:z], sr.gr[:ϵ], sr.gr[:def])
	elseif length(size(y)) == K2
		knots = (sr.gr[:b], sr.gr[:a], sr.gr[:z], sr.gr[:ϵ], sr.gr[:z], sr.gr[:ϵ])
	end

	itp_obj::itp_sr{length(size(y))} = interpolate(knots, y, Gridded(Linear()))

	return itp_obj
end