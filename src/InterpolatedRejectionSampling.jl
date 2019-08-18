#=
    InterpolatedRejectionSampling
    Copyright © 2019 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#

module InterpolatedRejectionSampling

using Interpolations
using Interpolations: Extrapolation
using NumericalIntegration
using StatsBase

export irsample!, irsample

include("envelope.jl")
include("irs_1dim.jl")

""" irsample(knots, probs, n)
inputs:
    knots -->  vector of knots
    probs -->  vector of probs
    n -->  number of samples to draw
"""
function irsample(knots::AbstractArray{T,D}, probs::AbstractArray{T,D}, n::Int
                 ) where {T<:Real,D}
    interp = LinearInterpolation(knots, probs, extrapolation_bc=Throw())
    return irsample(interp, n)
end

end
