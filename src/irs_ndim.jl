#=
    InterpolatedRejectionSampling
    Copyright © 2019 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#



""" function get_support(ks,ci)
input:
    ks -> tuple of knots (xknots, yknots, ...)
    ci -> CartesianIndex of bin
return:
    retval -> tuple of supports ((xmin,xspan), (ymin,yspan), ...)
"""
@inline function get_support(ks::NTuple{D,AbstractVector{T}},
                             ci::CartesianIndex{D}) where {D,T<:Real}
    return ntuple(i -> get_support(ks[i],ci.I[i]), Val(D))
end 



""" function get_sample(support)
input:
    support --> ((xmin,xspan),(ymin,yspan),...)
"""
@inline function get_sample(support::NTuple{D,NTuple{2,T}}) where {D,T<:Real}
    return [get_sample(s) for s in support]
end



""" function get_bin_knot_inds(ci::CartesianIndex{D})
get knot indices for ci bin
example: for CartesianIndex(2,3,5)
return (2:3, 3:4, 5:6)
"""
function get_bin_knot_inds(ci::CartesianIndex{D}) where {D}
    return ntuple(i -> ci.I[i]:ci.I[i]+1, Val(D))
end



""" function get_maxprob(ps,ci)
input:
    ps --> array of probs
    ci --> CartesianIndex of bin
"""
@inline function get_maxprob(ps::AbstractArray{T,D}, ci::CartesianIndex{D}
                            ) where {T<:Real,D}
    inds = get_bin_knot_inds(ci)
    return maximum(ps[CartesianIndices(inds)])
end 



""" function rsample(interp,cind)
input:
    interp --> interpolation object   
    ci --> CartesianIndex of bin to sample
"""
function rsample(interp::Extrapolation{T,D,ITPT,IT,ET}, cind::CartesianIndex{D}
                ) where {T<:Real,D,ITPT,IT,ET}

    support = get_support(get_knots(interp), cind)
    maxprob = get_maxprob(get_coefs(interp), cind)
    while true
        samp = get_sample(support)
        if rand()*maxprob ≤ interp(samp...)
            return samp
        end
    end
    error("unable to draw a sample after $maxruns runs")
end



""" function get_bin_ks_ps(interp, cind)
input:
    interp -->
    ci -->
"""
function get_bin_ks_ps(interp::Extrapolation{T,D,ITPT,IT,ET}, ci::CartesianIndex{D}
                      ) where {T<:Real,D,ITPT,IT,ET}
    inds = get_bin_knot_inds(ci)
    ks = ntuple(i -> get_knots(interp)[i][inds[i]], Val(D))
    ps = get_coefs(interp)[CartesianIndices(inds)]
    return ks, ps
end



""" function compute_pmf(interp)
input:
    interp
"""
function compute_pmf(interp::Extrapolation{T,D,ITPT,IT,ET}) where {T<:Real,D,ITPT,IT,ET}
    pmf = Array{T,D}(undef, size(get_coefs(interp)).-1 ...)
    for ci in CartesianIndices(pmf)
        bin_ks, bin_ps = get_bin_ks_ps(interp, ci)
        pmf[ci] = integrate(bin_ks, bin_ps)
    end
    return pmf./sum(pmf)
end



""" function choose_bins(interp, n)
"""
function choose_bins(interp::Extrapolation{T,D,ITPT,IT,ET}, n::Int
                    ) where {T<:Real,D,ITPT,IT,ET}
    pmf = compute_pmf(interp)
    return sample(CartesianIndices(pmf), Weights(vec(pmf)), n)
end

""" function choose_bin(interp)
"""
function choose_bin(interp::Extrapolation{T,D,ITPT,IT,ET}
                   ) where {T<:Real,D,ITPT,IT,ET}
    pmf = compute_pmf(interp)
    return sample(CartesianIndices(pmf), Weights(vec(pmf)))
end





""" function irsample(interp)
input:
    interp --> interpolation object
"""
function irsample(interp::Extrapolation{T,D,ITPT,IT,ET}
                 ) where {T<:Real,D,ITPT,IT,ET}

    bin = choose_bin(interp)
    return rsample(interp, bin)
end

""" function irsample(interp, n)
input:
    interp --> interpolation object
    n --> number of samples to draw
"""
function irsample(interp::Extrapolation{T,D,ITPT,IT,ET}, n::Int
                 ) where {T<:Real,D,ITPT,IT,ET}

    retval = Matrix{T}(undef,D,n)

    bins = choose_bins(interp, n)
    for (i,bin) in enumerate(bins)
        retval[:,i] = rsample(interp, bin)
    end
    return retval
end



""" irsample(knots, probs, n)
inputs:
    knots -->  vector of knots
    probs -->  vector of probs
    n -->  number of samples to draw
"""
function irsample(knots::NTuple{D,AbstractVector{T}},
                  probs::AbstractArray{T,D},
                  n::Int) where {T<:Real,D}
    interp = LinearInterpolation(knots, probs, extrapolation_bc=Throw())
    return irsample(interp, n)
end
