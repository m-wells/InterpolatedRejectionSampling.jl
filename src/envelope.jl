#=
    envelope
    Copyright © 2019 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#

""" function get_support(x)
input:
    x -> vector of knots
return:
    support -> (minimum(x), maximum(x) - minimum(x))
"""
@inline function get_support(knots ::AbstractVector{Float64}) ::NTuple{2,Float64}
    return (knots[1], knots[end]-knots[1])
end 

""" function get_support(X)
input:
    knots -> tuple of knots (xknots, yknots, ...)
return:
    retval -> tuple of supports ((xmin,xspan), (ymin,yspan), ...)
"""
@inline function get_support(knots ::NTuple{D,AbstractVector}) where D
    return ntuple(i -> get_support(knots[i]), Val(D))
end 

""" check_neighbors(x)
input:
    x -> boolean vector
output:
    boolean vector

returns true at indices where x[i-1] & x[i] & x[i+1] are true
"""
@inline function check_neighbors(x ::AbstractVector{Bool})
    retval = BitArray(x)
    if length(retval) > 1
        retval[1] &= x[2]
        retval[end] &= x[end-1]
        retval[2:end-1] .&= x[1:end-2] .& x[3:end]
    end
    return retval
end

""" check_neighbors(X)
input:
    X -> Boolean Array
output:
    Boolean Array

N-dimensional version of check_neighbors(x) 
"""
@inline function check_neighbors(X ::AbstractArray{Bool,D}) where D
    retval = BitArray(X)
    s = size(X)
    for i in 1:D
        for j in 1:s[i]
            _retval = selectdim(retval, i, j)
            _X = selectdim(X, i, j)
            _retval .&= check_neighbors(_X)
        end
    end
    return retval
end

""" get_support(knots, probs)
input:
    knots -> vectors specifing the knots
    probs -> vectors prob
output:
    min/span values

This routine accounts for range of zero prob and limits the support to exclude these regions.
Currently only outer zero regions can be excluded.
"""
@inline function get_support(knots ::AbstractVector,
                             probs ::AbstractVector)
    zmask = .!check_neighbors(probs .≤ 0)
    inds = findall(zmask)
    min_ind = minimum(inds)
    max_ind = maximum(inds)
    return get_support(knots[min_ind:max_ind])
end

""" get_support(knots, probs)
input:
    knots -> tuple of vectors specifing the knots
    probs -> array of prob
output:
    tuple of min/span values

This routine accounts for boxes of zero prob and limits the support to exclude these regions
Currently only outer zero regions can be excluded.
"""
@inline function get_support(knots ::NTuple{D,AbstractVector},
                     probs ::AbstractArray{Float64,D}) where D
    zmask = .!check_neighbors(iszero.(probs))
    cinds = findall(zmask)
    min_cind = minimum(cinds)
    max_cind = maximum(cinds)
    return ntuple(i -> get_support(knots[i][min_cind[i]:max_cind[i]]), Val(D))
end

""" get_bounding_box(pnt, knots)
input:
    pnt -> missing/float64 to surround if missing it is unbounded
    knots -> vector of knots
output:
    indices to access the the bounded range
"""
@inline function get_bounding_box(pnt ::Union{Missing,Float64},
                                  knots ::AbstractVector{Float64})
    left = Int(1)
    right = Int(length(knots))
    if !ismissing(pnt)
        left = Int(max(left, searchsortedlast(knots, pnt)))
        right = Int(min(right, searchsortedfirst(knots, pnt)))
    end
    return UnitRange(left,right)
end

""" get_bounding_box(pnt, knots)
inputs:
    pnt -> tuple of values (or missing) specify the point to bound
    knots -> tuple of knot vectors
output:
    tuple of indices that access the bounded box
"""
@inline function get_bounding_box(pnt ::NTuple{D,Union{Missing,Float64}},
                          knots ::NTuple{D,AbstractVector{Float64}}) where D
    return ntuple(i -> get_bounding_box(pnt[i], knots[i]), Val(D))
end

@inline function apply_bbox(x ::AbstractArray{T,N},
                            i ::NTuple{N,UnitRange{Int}}) where {T,N}
    return x[i...]
end

""" function get_support(bounding_box, knots, probs)
inputs:
    bounding_box -> indices of knots that specify the box
    knots -> tuple of knot vectors
    probs -> array of prob
output:
    support -> tuple of min/span
"""
@inline function get_support(bounding_box ::NTuple{D,UnitRange{Int}},
                             knots ::NTuple{D,AbstractVector},
                             probs ::AbstractArray{Float64,D}) where D
    _knots = ntuple(i -> knots[i][bounding_box[i]], Val(D))
    _probs = apply_bbox(probs, bounding_box)
    return get_support(_knots, _probs)
end

""" function get_envelope(pnt, knots, probs)
inputs:
    pnt -> tuple of floats/missing
    knots -> knots describing the grid of probs
    probs -> grid of probs
return: (support, pmax)
    support -> tuple of (xmin,xspan, ymin,yspan)
    pmax -> maximum prob accessible to pnt
"""
@inline function get_envelope(pnt ::NTuple{D,Union{Missing,Float64}},
                              knots ::NTuple{D,AbstractVector},
                              probs ::AbstractArray{Float64,D}) where D
    bounding_box = get_bounding_box(pnt, knots)
    support = get_support(bounding_box, knots, probs)
    pmax = maximum(apply_bbox(probs, bounding_box))
    return (support, pmax)
end

@inline function get_envelope(slices ::AbstractMatrix{Union{Missing,Float64}},
                              knots ::NTuple{D,AbstractVector},
                              probs ::AbstractArray{Float64,D}) where D
    n = last(size(slices))
    retval = [get_envelope(ntuple(i -> slices[i,j], Val(D)), knots, probs) for j = 1:n]
    return first.(retval), last.(retval)
end
