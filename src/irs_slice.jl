#=
    InterpolatedRejectionSampling
    Copyright © 2019 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#



function sliced_interp(slice::AbstractVector{Union{Missing,T}},
                       interp::Extrapolation{T,D,ITPT,IT,ET})

end


#function get_slice_ks_ps(s::AbstractVector{Union{Missing,T}},
#                         ks::NTuple{D,AbstractVector{T}},
#                         ps::AbstractArray{T,D}
#                        ) where {T<:Real,D}
#
#    inds = Vector{UnitRange{Int}}(undef,D)
#    for (i,si) in enumerate(s)
#        if ismissing(si)
#            inds[i] = 1:length(ks[i])
#        else
#            si < ks[i][1] && error("si < ks[i][1] encountered with si=",
#                                   si, " and ks[i][1]=",
#                                   ks[i][1])
#            si > ks[i][end] && error("si > ks[i][end] encountered with si=",
#                                     si, " and ks[i]=",
#                                     ks[i][end])
#            j = si == ks[i][1] ? 1 :
#                si == ks[i][end] ? length(ks[i]) - 1 :
#                searchsortedlast(ks[i], si)
#
#            inds[i] = j:j+1
#        end
#    end
#    s_ks = ntuple(i -> ks[i][inds[i]], Val(D))
#    s_ps = ps[CartesianIndices(tuple(inds...))]
#    return s_ks, s_ps
#end


""" function get_support_and_maxprob(s, ks, ps, ci)
input:
    s -> slice vector
    ks -> tuple of knot vectors
    ci -> CartesianIndex of bin
"""
function get_support(s::AbstractVector{Union{Missing,T}},
                     ks::NTuple{D,AbstractVector{T}},
                     ci::CartesianIndex{D}) where {D,T<:Real}
    return ntuple(i -> ismissing(s[i]) ? missing : get_support(ks[i],ci.I[i]), Val(D))
end

""" function get_support_and_maxprob(s,ks,ps,ci)
input:
    s -> slice vector
    ks -> tuple of knot vectors
    ps -> probability array
    ci -> CartesianIndex of bin
"""
function get_support_and_maxprob(s::AbstractVector{Union{Missing,T}},
                                 interp::Extrapolation,
                                 ci::CartesianIndex{D}) where {D,T<:Real}
    inds = get_bin_inds(ci)
    probs = ps[CartesianIndices(inds)]
    return get_support(ks,ci), maximum(probs)
end 



""" function get_bin_ks_ps(slice, interp)
input:
    slice -->
    interp -->
    ci -->
"""
function get_bin_ks_ps(slice::AbstractVector{Union{Missing,T}},
                       interp::Extrapolation{T,D,ITPT,IT,ET},
                       ci::CartesianIndex{D}
                      ) where {T<:Real,D,ITPT,IT,ET}
    inds = get_bin_knot_inds(ci)
    ks = ntuple(i -> interp.itp.knots[i][inds[i]], Val(D))
    ps = interp.itp.coefs[CartesianIndices(inds)]
    return ks, ps
end




""" function compute_pmf(slice,interp)
input:
    slice
    interp
"""
function compute_pmf(slice::AbstractVector{Union{Missing,T}},
                     interp::Extrapolation{T,D,ITPT,IT,ET}) where {T<:Real,D,ITPT,IT,ET}
    mask = ismissing.(slice)
    dim_size = size(interp.itp.coefs)[mask] .-1
    dim = length(dim_size)

    pmf = Array{T,dim}(undef, dim_size...)
    for ci in CartesianIndices(pmf)
        bin_ks, bin_ps = get_bin_ks_ps(slice, interp, ci)
        pmf[ci] = integrate(bin_ks, bin_ps)
    end
    return pmf./sum(pmf)
end



""" function choose_bin(slice, interp)
"""
function choose_bin(slice::AbstractVector{Union{Missing,T}},
                    interp::Extrapolation{T,D,ITPT,IT,ET}
                   ) where {T<:Real,D,ITPT,IT,ET}
    pmf = compute_pmf(slice, interp)
    return sample(CartesianIndices(pmf), Weights(vec(pmf)))
end



""" irsample(slice, knots, probs)
inputs:
    slices -->  vector of dimensionality D with missing values to be drawn
    interp --> interpolation object
"""
function irsample!(slices::Matrix{Union{Missing,T}},
                   interp::Extrapolation{T,D,ITPT,IT,ET}
                 ) where {T<:Real,D,ITPT,IT,ET}

    # find columns in slices which contain any missing values
    slice_mask = [any(ismissing.(s)) for s in eachcol(slices)]
    slice_inds = findall(slice_mask)

    for i in slice_inds
        slice = selectdim(slices, 2, i)
        s_ks, s_ps = get_slice_ks_ps(slice, knots, probs)
        bin = choose_bin(s_ks, s_ps)

        slices[:,i] = rsample(slice, interp, bin)
    end
end
