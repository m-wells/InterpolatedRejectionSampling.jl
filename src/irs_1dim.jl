#=
    InterpolatedRejectionSampling
    Copyright © 2019 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#



""" function get_support(x,i)
input:
    x --> vector of knots
    i --> bin number
"""
@inline function get_support(x::AbstractVector{T}, i::Int) where T<:Real
    return (x[i], x[i+1] - x[i])
end 

@inline function get_support(x::AbstractVector{T}, ci::CartesianIndex{1}) where T<:Real
    i = first(ci.I)
    return (x[ci], x[i+1] - x[i])
end 



""" function get_sample(support)
input:
    support --> the minimum and span of the region being sampled
"""
@inline function get_sample(support::NTuple{2,T}) where T<:Real
    return support[1] + rand()*support[2]
end



""" function get_maxprob(p,i)
input:
    p --> vector of prob
    i --> bin number
"""
@inline function get_maxprob(p::AbstractVector{T}, i::Int) where T<:Real
    return max(p[i], p[i+1])
end 



""" function rsample(interp,i)
input:
    interp --> interpolation object   
    i --> bin to sample
"""
function rsample(interp::Extrapolation{T,1,ITPT,IT,ET}, i::Int
                ) where {T<:Real,ITPT,IT,ET}
    support = get_support(get_knots(interp), i)
    maxprob = get_maxprob(get_coefs(interp), i)
    while true
        samp = get_sample(support)
        if rand()*maxprob ≤ interp(samp)
            return samp
        end
    end
    error("unable to draw a sample after $maxruns runs")
end



""" function compute_pmf(interp)
Computes the probability mass function for each bin
input:
    interp --> interpolation object
"""
function compute_pmf(interp::Extrapolation{T,1,ITPT,IT,ET}
                    ) where {T<:Real,ITPT<:GriddedInterpolation,IT,ET}
    # trapezoidal rule
    pmf = diff(get_knots(interp)).*midpoints(get_coefs(interp)).*(1//2)
    return pmf./sum(pmf)
end

function compute_pmf(interp::Extrapolation{T,1,ITPT,IT,ET}
                    ) where {T<:Real,ITPT<:ScaledInterpolation,IT,ET}
    # trapezoidal rule
    pmf = step(get_knots(interp)).*midpoints(get_coefs(interp)).*(1//2)
    return pmf./sum(pmf)
end





""" function choose_bins(interp, n)
input:
    interp --> interpolation object
    n --> number of bins to choose
"""
function choose_bins(interp::Extrapolation{T,1,ITPT,IT,ET}, n::Int) where {T<:Real,ITPT,IT,ET}
    pmf = compute_pmf(interp)
    return sample(eachindex(pmf), Weights(pmf), n)
end



""" function irsample(interp, n)
input:
    interp --> interpolation object
    n --> number of samples to draw
"""
function irsample(interp::Extrapolation{T,1,ITPT,IT,ET}, n::Int
                 ) where {T<:Real,ITPT,IT,ET}

    retval = Vector{T}(undef,n)

    bins = choose_bins(interp, n)
    for (i,bin) in enumerate(bins)
        retval[i] = rsample(interp, bin)
    end
    return retval
end

""" irsample(knots, probs, n)
inputs:
    knots -->  vector of knots
    probs -->  vector of probs
    n -->  number of samples to draw
"""
function irsample(knots::AbstractVector{T},
                  probs::AbstractVector{T},
                  n::Int) where {T<:Real}
    interp = LinearInterpolation(knots, probs, extrapolation_bc=Throw())
    return irsample(interp, n)
end
