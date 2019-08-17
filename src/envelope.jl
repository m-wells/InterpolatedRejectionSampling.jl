#=
    envelope
    Copyright © 2019 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#

"""
compute_pmf

Computes the probability mass function for each bin
"""
@inline function compute_pmf(knots::AbstractVector{T},
                             probs::AbstractVector{T}
                            ) where T<:Real
    # trapezoidal rule
    pmf = diff(knots).*midpoints(probs).*(1//2)
    return pmf./sum(pmf)
end

"""
get_bins

Returns bins sampled using the computed pmf
"""
@inline function get_bins(k::AbstractVector{T}, p::AbstractVector{T}, n::Int) where T<:Real
    # compute the probability mass function for each bin
    pmf = compute_pmf(k, p)
    # select bins 
    w = Weights(pmf)
    return sample(eachindex(pmf), w, n)
end

""" function get_bin_inds(ci::CartesianIndex{D})
for CartesianIndex(2,3,5)
return (2:3, 3:4, 5:6)
"""
function get_bin_inds(ci::CartesianIndex{D}) where {D}
    return ntuple(i -> ci.I[i]:ci.I[i]+1, D)
end

function compute_pmf(knots::NTuple{D,AbstractVector{T}},
                     probs::AbstractArray{T,D},
                    ) where {D,T<:Real}
    pmf = Array{T}(undef,size(probs) .-1 ...)
    for ci in CartesianIndices(pmf)
        inds = get_bin_inds(ci)
        cknots = ntuple(i -> knots[i][inds[i]], D)
        cprobs = probs[CartesianIndices(inds)]
        pmf[ci] = integrate(cknots, cprobs)
    end
    return pmf./sum(pmf)
end

function get_bins(knots::NTuple{D,AbstractVector{T}},
                  probs::AbstractArray{T,D},
                  n::Int) where {D,T<:Real}
    # compute the probability mass function for each bin
    pmf = compute_pmf(knots, probs)
    # select bins 
    w = Weights(vec(pmf))
    return sample(CartesianIndices(pmf),w,n)
end
function get_bin(knots::NTuple{D,AbstractVector{T}},
                 probs::AbstractArray{T,D},
                ) where {D,T<:Real}
    return get_bins(knots, probs, 1)[1]
end

""" function get_support(x)
input:
    x -> vector of knots
    i -> bin number
return:
    support and maxprob
"""
@inline function get_support(x::AbstractVector{T}, i::Int) where T<:Real
    return (x[i], x[i+1] - x[i])
end 

""" function get_support_and_maxprob(x,p,i)
input:
    x -> vector of knots
    p -> probability at knots
    i -> bin number
"""
function get_support_and_maxprob(x::AbstractVector{T},
                                 p::AbstractVector{T},
                                 i::Int) where T<:Real
    return get_support(x,i), max(p[i],p[i+1])
end

""" function get_support(ks,ci)
input:
    ks -> tuple of knots (xknots, yknots, ...)
    ci -> CartesianIndex of bin
return:
    retval -> tuple of supports ((xmin,xspan), (ymin,yspan), ...)
"""
@inline function get_support(ks::NTuple{D,AbstractVector{T}},
                             ci::CartesianIndex{D}) where {D,T<:Real}
    return ntuple(i -> get_support(ks[i],ci.I[i]), D)
end 

""" function get_support_and_maxprob(ks,ps,ci)
input:
    ks -> tuple of knot vectors
    ps -> probability array
    ci -> CartesianIndex of bin
"""
@inline function get_support_and_maxprob(ks::NTuple{D,AbstractVector{T}},
                                         ps::Array{T,D},
                                         ci::CartesianIndex{D}) where {D,T<:Real}
    inds = get_bin_inds(ci)
    probs = ps[CartesianIndices(inds)]
    return get_support(ks,ci), maximum(probs)
end 
