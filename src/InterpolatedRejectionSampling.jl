#=
    InterpolatedRejectionSampling
    Copyright © 2019 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#

module InterpolatedRejectionSampling

using Interpolations

export irsample!, irsample

include("envelope.jl")

#const maxruns = 1_000_000_000
#const maxruns = 10

#@inline function get_sample!(cache ::Vector{Float64},
#                              pnt ::Union{Missing,Float64},
#                              support ::NTuple{2,Float64})
#    if ismissing(pnt)
#        return support[1] + popfirst!(cache)*support[2]
#    else
#        return pnt
#    end
#end

@inline function get_sample(support ::NTuple{2,Float64}) where D
    return support[1] + rand()*support[2]
end

@inline function get_sample(support ::NTuple{D,NTuple{2,Float64}}) where D
    return Float64[ntuple(i -> get_sample(support[i]), Val(D))...]
end

@inline function get_sample(slice ::AbstractVector{Union{Missing,Float64}},
                            support ::NTuple{D,NTuple{2,Float64}}) where D
    
    return Float64[ntuple(i -> ismissing(slice[i]) ? get_sample(support[i]) : slice[i], Val(D))...]
end

@inline function rsample(interp ::Interpolations.Extrapolation,
                         support ::NTuple{2,Float64},
                         maxprob ::Float64) where D
    while true
        samp = get_sample(support)
        if rand()*maxprob ≤ interp(samp)
            return samp
        end
    end
    error("unable to draw a sample after $maxruns runs")
end

@inline function rsample(interp ::Interpolations.Extrapolation,
                         support ::NTuple{D,NTuple{2,Float64}},
                         maxprob ::Float64) where D
    while true
        samp = get_sample(support)
        if rand()*maxprob ≤ interp(samp...)
            return samp
        end
    end
    error("unable to draw a sample after $maxruns runs")
end

@inline function rsample(interp ::Interpolations.Extrapolation,
                         slice ::AbstractVector{Union{Missing,Float64}},
                         support ::NTuple{D,NTuple{2,Float64}},
                         maxprob ::Float64) where D

    while true
        samp = get_sample(slice, support)
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
@inline function irsample(knots ::AbstractVector{Float64},
                          probs ::AbstractVector{Float64},
                          n ::Int)
    retval = Vector{Float64}(undef, n)
    interp = LinearInterpolation(knots, probs, extrapolation_bc=Throw())

    support = get_support(knots)
    maxprob = maximum(probs)

    for i in 1:n
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
@inline function irsample(knots ::NTuple{D,AbstractVector{Float64}},
                          probs ::AbstractArray{Float64,D},
                          n ::Int) where D
    retval = Matrix{Float64}(undef, D, n)
    interp = LinearInterpolation(knots, probs, extrapolation_bc=Throw())

    support = get_support(knots)
    maxprob = maximum(probs)
    for i in 1:n
        retval[:,i] = rsample(interp, support, maxprob)
    end
    return retval
end

""" irsample!(slices, knots, probs)
inputs:
    slices ⟶  D by N matrix where D is the dimensionality and N is the number of points, each column is a point (xi,yi,...), missing values will be drawn
    knots ⟶  tuple of knots (xknots, yknots, ...)
    probs ⟶  D-dimension grid of probability
"""
@inline function irsample!(slices ::Matrix{Union{Missing,Float64}},
                           knots ::NTuple{D,AbstractVector{Float64}},
                           probs ::AbstractArray{Float64,D}) where D
    interp = LinearInterpolation(knots, probs, extrapolation_bc=Throw())
    nslices = last(size(slices))
    mask = ismissing.(slices)
    slice_inds = findall([any(mask[:,i]) for i = 1:nslices])

    support = get_support(knots)
    maxprob = maximum(probs)
    for (i,ind) in enumerate(slice_inds)
        slice = selectdim(slices, 2, ind)
        slice .= rsample(interp, vec(slice), support, maxprob)
    end
end

end

#""" irsample!(slices, knots, probs)
#inputs:
#    slices ⟶  D by N matrix where D is the dimensionality and N is the number of points, each column is a point (xi,yi,...), missing values will be drawn
#    knots ⟶  tuple of knots (xknots, yknots, ...)
#    probs ⟶  D-dimension grid of probability
#"""
#@inline function irsample2!(slices ::Matrix{Union{Missing,Float64}},
#                           knots ::NTuple{D,AbstractVector{Float64}},
#                           probs ::AbstractArray{Float64,D};
#                           verbose ::Bool = false) where D
#
#    interp = LinearInterpolation(knots, probs, extrapolation_bc=Throw())
#
#    nslices = last(size(slices))
#    mask = ismissing.(slices)
#    slice_inds = findall([any(mask[:,i]) for i = 1:nslices])
#    nmissing = [count(mask[:,i]) for i = 1:nslices]
#
#    support, pmaxes = get_envelope(slices, knots, probs)
#
#    for run = 0:maxruns
#        #@info "slice_inds" slice_inds
#        _slices = selectdim(slices, 2, slice_inds)
#        _support = view(support, slice_inds)
#        _pmaxes = view(pmaxes, slice_inds)
#        if length(slice_inds) == 0
#            return
#        end
#        #if verbose & (mod(run,10) == 0)
#        if verbose
#            @info """
#            running "irsample!" on run $run with $(length(slice_inds)) to draw
#            """
#        end
#
#        #@info "generating cache " length(slice_inds) nmissing[slice_inds]
#        cache = rand(length(slice_inds) + sum(nmissing[slice_inds]))
#
#        samples = Matrix{Float64}(undef,D,length(slice_inds))
#        for ind in eachindex(_slices)
#            samples[ind] = get_sample!(cache, _slices, ind, _support)
#        end
#
#        variate = cache.*_pmaxes
#
#        inds = findall([variate[ind] ≤ interp(samples[:,ind]...) for ind = 1:length(slice_inds)])
#        # update slices
#        _slices[:,inds] .= samples[:,inds]
#        deleteat!(slice_inds, inds)
#
#        #slice_inds = findall([any(ismissing.(_slices[:,i])) for i = 1:size(_slices)[2]])
#        #_slices = selectdim(slices, 2, slice_inds)
#        #_nslices = length(slice_inds)
#    end
#    nmissing = sum(count.([ismissing.(s) for s in slices]))
#    if nmissing != 0
#        error("failed to draw all samples after maximum number of runs") 
#    end
#end





#""" irsample!(slices, knots, probs)
#inputs:
#    slices ⟶  D by N matrix where D is the dimensionality and N is the number of points, each column is a point (xi,yi,...), missing values will be drawn
#    knots ⟶  tuple of knots (xknots, yknots, ...)
#    probs ⟶  D-dimension grid of probability
#"""
#function irsample!(slices ::AbstractMatrix{Union{Missing,Float64}},
#                   knots ::NTuple{D,AbstractVector{Float64}},
#                   probs ::AbstractArray{Float64,D};
#                   verbose ::Bool = false) where D
#
#    interp = LinearInterpolation(knots, probs, extrapolation_bc=Throw())
#    support, pmaxes = get_envelope(slices, knots, probs)
#    pmaxes .+= eps.(pmaxes)
#
#    _, nc = size(slices)
#
#    for run = 1:maxruns
#        mask = BitArray(ismissing.(slices))
#        nmissing = count(mask)
#        if verbose & mod(run,10) == 0
#            @show run
#            @show nmissing
#        end
#
#        if nmissing == 0
#            break
#        end
#
#        slice_inds = findall([any(mask[:,i]) for i = 1:nc])
#        nslices = length(slice_inds)
#
#        cache = rand(nslices + nmissing)
#
#        _slices = selectdim(slices,2,slice_inds)
#
#        _support = view(support,slice_inds)
#        _pmaxes = view(pmaxes,slice_inds)
#        samples = collect(_slices)
#
#        cinds = findall(ismissing.(samples))
#        for (i,cind) in enumerate(cinds)
#            dim = cind[1]
#            ind = cind[2]
#            s_min = _support[ind][dim][1]
#            s_span = _support[ind][dim][2]
#            samples[cind] = s_min + pop!(cache)*s_span
#        end
#
#        variate = cache.*_pmaxes
#
#        inds = findall([variate[i] ≤ interp(samples[:,i]...) for i in eachindex(variate)])
#        # update slices
#        _slices[:,inds] = samples[:,inds]
#    end
#
#    mask = BitArray(ismissing.(slices))
#    nmissing = count(mask)
#    if nmissing != 0
#        error("failed to draw all samples after maximum number of runs") 
#    end
#end


