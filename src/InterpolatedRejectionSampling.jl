#=
    InterpolatedRejectionSampling
    Copyright © 2019 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#

module InterpolatedRejectionSampling

using Interpolations

export irsample!, irsample

"""
univariate version
`support` is of the form `(x1,xspan)` where `xspan` is `x2 - x1`.
"""
function get_support(knots :: AbstractVector)
    return (knots[1],knots[end] - knots[1])
end 

""" function get_support(knots)
input:
    knots ⟶  tuple of knots (xknots, yknots, ...)
return:
    retval ⟶  tuple of supports ((xmin,xspan), (ymin,yspan), ...)
"""
function get_support(knots :: NTuple{D,AbstractVector}) where D
    return ntuple(i -> (get_support(knots[i])), Val(D))
end 

""" function get_pmax(pnt, knots, probs)
input:
    pnt ⟶  tuple of coordinates (x,y,...) with missing values acting as Colon()
    knots ⟶  tuple of knots (xknots, yknots, ...)
    probs ⟶  D-dimension grid of probability
return:
    pmax ⟶  largest value of probs that pnt could take
"""
function get_pmax(pnt ::NTuple{D,Union{Missing,Float64}},
                  knots ::NTuple{D,AbstractVector{Float64}},
                  probs ::AbstractArray{Float64,D}) where D
    pnts = ntuple(i -> ifelse(ismissing(pnt[i]), 1:length(knots[i]),
                              max(1,searchsortedlast(knots[i],pnt[i])):min(length(knots[i]),searchsortedfirst(knots[i],pnt[i]))),
                  Val(D))
    cartinds = CartesianIndices(pnts)
    return maximum(probs[cartinds])
end

""" irsample!(slices, knots, probs)
inputs:
    slices ⟶  D by N matrix where D is the dimensionality and N is the number of points, each column is a point (xi,yi,...), missing values will be drawn
    knots ⟶  tuple of knots (xknots, yknots, ...)
    probs ⟶  D-dimension grid of probability
"""
function irsample!(slices ::AbstractMatrix{Union{Missing,Float64}},
                   knots ::NTuple{D,AbstractVector{Float64}},
                   probs ::AbstractArray{Float64,D}) where D

    interp = LinearInterpolation(knots, probs, extrapolation_bc=Throw())
    support = get_support(knots)

    _, nc = size(slices)

    for run = 1:100_000
        mask = BitArray(ismissing.(slices))
        nmissing = count(mask)

        if nmissing == 0
            break
        end

        inds = findall([any(mask[:,i]) for i = 1:nc])
        nslices = length(inds)

        cache = rand(nslices + nmissing)

        _slices = selectdim(slices,2,inds)
        samples = collect(_slices)

        pmaxes = [get_pmax(ntuple(i->samples[i,j], Val(D)), knots, probs) for j = 1:nslices]
        pmaxes .+= eps.(pmaxes)

        cinds = findall(ismissing.(samples))
        for cind in cinds
            samples[cind] = support[cind[1]][1] + pop!(cache)*support[cind[1]][2]
        end

        variate = cache.*pmaxes

        inds = findall([variate[i] ≤ interp(samples[:,i]...) for i in eachindex(variate)])
        # update slices
        _slices[:,inds] = samples[:,inds]
    end

    mask = BitArray(ismissing.(slices))
    nmissing = count(mask)
    if nmissing != 0
        error("failed to draw all samples after maximum number of runs") 
    end
end

""" irsample(knots, probs, n)
inputs:
    knots ⟶  tuple of knots (xknots, yknots, ...)
    probs ⟶  D-dimension grid of probability
    n ⟶  number of samples to draw
"""
function irsample(knots ::NTuple{D,AbstractVector{Float64}},
                  probs ::AbstractArray{Float64,D},
                  n ::Int) where D

    retval = Matrix{Union{Missing,Float64}}(missing, D, n)
    irsample!(retval, knots, probs)
    return convert(Matrix{Float64}, retval)
end

""" irsample(knots, probs, n)
inputs:
    knots ⟶  vector of knots
    probs ⟶  vector of probs
    n ⟶  number of samples to draw
"""
function irsample(knots ::AbstractVector{Float64},
                  probs ::AbstractVector{Float64},
                  n ::Int)

    retval = Matrix{Union{Missing,Float64}}(missing, 1, n)
    irsample!(retval, (knots,), probs)
    return convert(Vector{Float64}, vec(retval))
end

end
