#=
    InterpolatedRejectionSampling
    Copyright © 2019 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#

module InterpolatedRejectionSampling

using Interpolations
using NumericalIntegration
using StatsBase

export irsample!, irsample

include("envelope.jl")

@inline function get_sample(support::NTuple{2,T}) where T<:Real
    return support[1] + rand()*support[2]
end

@inline function get_sample(support::NTuple{D,NTuple{2,T}}) where {D,T<:Real}
    return [get_sample(s) for s in support]
end

#@inline function get_sample(slice ::AbstractVector{Union{Missing,Float64}},
#                            support ::NTuple{D,NTuple{2,Float64}}) where D
#    return Float64[ntuple(i -> ismissing(slice[i]) ? get_sample(support[i]) : slice[i], Val(D))...]
#end

@inline function rsample(interp ::Interpolations.Extrapolation,
                         support ::NTuple{2,T},
                         maxprob ::T) where {D,T<:Real}
    while true
        samp = get_sample(support)::T
        if rand()*maxprob ≤ interp(samp)
            return samp
        end
    end
    error("unable to draw a sample after $maxruns runs")
end

@inline function rsample(interp ::Interpolations.Extrapolation,
                         support ::NTuple{D,NTuple{2,T}},
                         maxprob ::T) where {D,T<:Real}
    while true
        samp = get_sample(support)::Vector{T}
        if rand()*maxprob ≤ interp(samp...)
            return samp
        end
    end
    error("unable to draw a sample after $maxruns runs")
end

""" irsample(knots, probs, n)
inputs:
    knots ⟶  vector of knots
    probs ⟶  vector of probs
    n ⟶  number of samples to draw
"""
@inline function irsample(knots ::AbstractVector{T},
                          probs ::AbstractVector{T},
                          n ::Int) where T<:Real

    retval = Vector{T}(undef, n)
    interp = LinearInterpolation(knots, probs, extrapolation_bc=Throw())
    bins = get_bins(knots, probs, n)

    for (i,bin) in enumerate(bins)
        support, maxprob = get_support_and_maxprob(knots, probs, bin)
        retval[i] = rsample(interp, support, maxprob)
    end
    return retval
end

""" irsample(knots, probs, n)
inputs:
    knots ⟶  tuple of knots (xknots, yknots, ...)
    probs ⟶  D-dimension grid of probability
    n ⟶  number of samples to draw
"""
@inline function irsample(knots ::NTuple{D,AbstractVector{T}},
                          probs ::Array{T,D},
                          n ::Int) where {D,T<:Real}
    retval = Matrix{T}(undef, D, n)
    interp = LinearInterpolation(knots, probs, extrapolation_bc=Throw())
    bins = get_bins(knots, probs, n)
    for (i,bin) in enumerate(bins)
        support, maxprob = get_support_and_maxprob(knots, probs, bin)
        retval[:,i] = rsample(interp, support, maxprob)
    end
    return retval
end

function get_slice_knots(slice::AbstractVector{Union{Missing,T}},
                         knots::NTuple{D,AbstractVector{T}}
                        ) where {D,T<:Real}
    return ntuple(i -> ismissing(slice[i]) ? knots[i] : [slice[i]], D)
end

function get_slice_probs(s_knots::NTuple{D,AbstractVector{T}},
                         interp::Interpolations.Extrapolation
                        ) where {D,T<:Real}
    s_probs = Array{T}(undef, length.(s_knots)...)
    for ci in CartesianIndices(s_probs)
        #knot = ntuple(i -> s_knots[i][ci.I[i]], D)::NTuple{D,T}
        knot = [s_knots[i][ci.I[i]] for i in 1:D]
        s_probs[ci] = interp(knot...)::T
    end
    return s_probs
end

""" irsample(slice, knots, probs)
inputs:
    slice ⟶  vector of dimensionality D with missing values to be drawn
    knots ⟶  tuple of knots (xknots, yknots, ...)
    probs ⟶  D-dimension grid of probability
"""
@inline function irsample!(slices::Matrix{Union{Missing,T}},
                           knots::NTuple{D,AbstractVector{T}},
                           probs::Array{T,D}) where {D,T<:Real}
    interp = LinearInterpolation(knots, probs, extrapolation_bc=Throw())
    nslices = last(size(slices))

    slice_mask = [any(ismissing.(s)) for s in eachcol(slices)]
    slice_inds = findall(slice_mask)

    for i in slice_inds
        slice = selectdim(slices, 2, i)
        s_knots = get_slice_knots(slice, knots)
        s_probs = get_slice_probs(s_knots, interp)::Array{T,D}
        bin = get_bin(s_knots, s_probs)::CartesianIndex{D}

        support, maxprob = get_support_and_maxprob(s_knots, s_probs, bin)
        slices[:,i] = rsample(interp, support, maxprob)::Vector{T}
    end
end

end
