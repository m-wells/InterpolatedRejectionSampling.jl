""" function sliced_interp(slice, interp)
input:
    slice --> Slice "object"
    interp --> parent interpolation
"""
function sliced_interp(slice::Slice{T,D},
                       interp::Extrapolation{T,D,ITPT,IT,ET},
                      ) where {T<:Real,D,ITPT,IT,ET}
    d = length(size(slice))
    coefs = Array{T,d}(undef,size(slice)...)
    for (i,knot) in enumerate(slice)
        coefs[i] = interp(knot...)
    end
    knots = slice.knots[slice.mask]
    return LinearInterpolation(knots, coefs)
end





""" function common_mask_check(slices)
"""
function common_mask_check(slices::Matrix{Union{Missing,T}}) where {T<:Real}
    nrows = first(size(slices))
    rowcheck = Vector{Bool}(undef,nrows)

    #for (i,x) in enumerate(eachrow(slices))
    for i in 1:nrows
        # all or nothing check
        rowcheck[i] = all(ismissing.(slices[i,:])) || !any(ismissing.(slices[i,:]))
    end

    return all(rowcheck)
end



""" function irsample!(slices, interp)
"""
function irsample!(slices::Matrix{Union{Missing,T}},
                   interp::Extrapolation{T,D,ITPT,IT,ET}
                  ) where {T<:Real,D,ITPT,IT,ET}
    #for (i,s) in enumerate(eachcol(slices))
    ncols = last(size(slices))
    for i in 1:ncols
        s = view(slices, :, i)
        slice = Slice(s,interp)
        sinterp = sliced_interp(slice,interp)
        slices[slice.mask, i] .= irsample(sinterp)
    end
end



""" function irsample!(slices, interp, mask)
"""
function irsample!(slices::Matrix{Union{Missing,T}},
                   interp::Extrapolation{T,D,ITPT,IT,ET},
                   mask::BitVector
                  ) where {T<:Real,D,ITPT,IT,ET}

    cinds = CartesianIndices(eachindex.(get_knots(interp)[mask]))
    #for (i,s) in enumerate(eachcol(slices))
    ncols = last(size(slices))
    for i in 1:ncols
        s = view(slices, :, i)
        slice = Slice(s,interp,mask,cinds)
        sinterp = sliced_interp(slice,interp)
        temp = irsample(sinterp)
        slices[mask, i] .= irsample(sinterp)
    end
end



""" function irsample!(slices, knots, probs)
"""
function irsample!(slices::Matrix{Union{Missing,T}},
                   knots::NTuple{D,AbstractVector{T}},
                   probs::AbstractArray{T,D}
                  ) where {T<:Real,D}
    interp = LinearInterpolation(knots, probs)
    # check for a common masking
    if common_mask_check(slices)
        mask = ismissing.(slices[:,1])
        irsample!(slices, interp, mask)
    else
        irsample!(slices, interp)
    end
end
