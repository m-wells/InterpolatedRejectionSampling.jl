#=
    InterpolatedRejectionSampling
    Copyright © 2019 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#

module InterpolatedRejectionSampling

using Interpolations
using Interpolations: Extrapolation, GriddedInterpolation
using StatsBase

include("NumericalIntegration.jl")
using InterpolatedRejectionSampling.NumericalIntegration

export irsample!, irsample



struct Slice{T,D}
    slice::Vector{Union{Missing,T}}
    knots::NTuple{D,AbstractVector{T}}
    cinds::CartesianIndices
    mask::BitVector

    function Slice(slice::AbstractVector{Union{Missing,T}},
                   knots::NTuple{D,AbstractVector{T}},
                   mask  = ismissing.(slice)::BitVector,
                   cinds = CartesianIndices(eachindex.(knots[mask]))
                  ) where {T<:Real,D}
        return new{T,D}(slice,knots,cinds,mask)
    end

    function Slice(slice::AbstractVector{Union{Missing,T}},
                   interp::Extrapolation{T,D,ITPT,IT,ET},
                   mask  = ismissing.(slice)::BitVector,
                   cinds = CartesianIndices(eachindex.(get_knots(interp)[mask]))
                  ) where {T<:Real,D,ITPT,IT,ET}
        return new{T,D}(slice,get_knots(interp),cinds,mask)
    end
end



function Base.iterate(S::Slice{T,D}, state=1) where {T,D}
    state > length(S.cinds) && return nothing

    cinds = Vector{Int}(undef, D)
    cinds[S.mask] .= S.cinds[state].I
    return (ntuple(i -> S.mask[i] ? S.knots[i][cinds[i]] : S.slice[i], Val(D)), state+1)
end



function Base.size(S::Slice{T,D}) where {T,D}
    return length.(S.knots[S.mask])
end



@inline get_knots(interp::Extrapolation{T,1,ITPT,IT,ET}
                 ) where {T,D,ITPT<:GriddedInterpolation,IT,ET} = first(interp.itp.knots)
@inline get_knots(interp::Extrapolation{T,1,ITPT,IT,ET}
                 ) where {T,D,ITPT<:ScaledInterpolation,IT,ET} = first(interp.itp.ranges)
@inline get_knots(interp::Extrapolation{T,D,ITPT,IT,ET}
                 ) where {T,D,ITPT<:GriddedInterpolation,IT,ET} = interp.itp.knots
@inline get_knots(interp::Extrapolation{T,D,ITPT,IT,ET}
                 ) where {T,D,ITPT<:ScaledInterpolation,IT,ET} = interp.itp.ranges

@inline get_coefs(interp::Extrapolation{T,D,ITPT,IT,ET}
                 ) where {T,D,ITPT<:GriddedInterpolation,IT,ET} = interp.itp.coefs
@inline get_coefs(interp::Extrapolation{T,D,ITPT,IT,ET}
                 ) where {T,D,ITPT<:ScaledInterpolation,IT,ET} = interp.itp.itp.coefs


include("irs_1dim.jl")
include("irs_ndim.jl")
include("irs_slice.jl")

end
