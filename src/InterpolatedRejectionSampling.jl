#=
    InterpolatedRejectionSampling
    Copyright © 2019 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#

module InterpolatedRejectionSampling
#---------------------------------------------------------------------------------------------------
using Interpolations
####################################################################################################
#
####################################################################################################
function get_variate(support :: NTuple{N,NTuple{2,T}}) where {T,N}
    return ntuple( i -> (support[i][1] + rand(T)*support[i][2]), Val(N))
end
####################################################################################################
#
####################################################################################################
function rsample( sitp :: ScaledInterpolation{T,N,ITPT,IT,RT}
                , support :: NTuple{N,NTuple{2,T}}
                , pmax :: T
                ) where {T,N,ITPT,IT,RT}
    # draw a sample
    samp = get_variate(support)
    # determine a rejection value
    prob = rand()*pmax

    # if rejection value is greater than the interpolated value of the function at sample
    while prob > sitp(samp...)
        # redraw
        samp = get_variate(support)
        prob = rand()*pmax
    end
    return samp
end
####################################################################################################
#
####################################################################################################
function rejection_sampling( n :: Integer
                           , pdf   :: AbstractArray{T,N}
                           , knots :: Vararg{AbstractRange,N}
                           ) where {T,N}
                           #) :: Vector{NTuple{N,Float64}} where {T,N}
    # interpolate the pdf with cubic spline with linear edges with knots on cell edges
    itp = interpolate(pdf,BSpline(Cubic(Line(OnGrid()))))
    sitp = scale(itp, knots...)
    support = ntuple( i -> ( sitp.ranges[i][1]
                            , sitp.ranges[i][end] - sitp.ranges[i][1]
                           )
                    , Val(N)
                    )

    samp = Vector{NTuple{N,T}}(undef,n)
    pmax = 1.05*maximum(pdf)
    for i in eachindex(samp)
        samp[i] = rsample(sitp,support,pmax)
    end

    return samp
end
####################################################################################################
#
####################################################################################################
#---------------------------------------------------------------------------------------------------
end
