#=
    InterpolatedRejectionSampling
    Copyright © 2019 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#

module InterpolatedRejectionSampling

using Interpolations

export irsample, irsample!

const PROB_PAD = 3//2

abstract type InterpolationMethod end
struct Linear <: InterpolationMethod end
struct CubicSpline <: InterpolationMethod end

"""
n-dimensional version
`support` is of the form `((x1,xspan),(y1,yspan),...)` where `xspan` is `x2 - x1`.
"""
function get_support(knots :: NTuple{N,AbstractVector}) where N
    return ntuple( i -> ( knots[i][1]
                        , knots[i][end] - knots[i][1]
                        )
                 , Val(N)
                 )
end 

"""
univariate version
`support` is of the form `(x1,xspan)` where `xspan` is `x2 - x1`.
"""
function get_support(knots :: AbstractVector)
    return (knots[1],knots[end] - knots[1])
end 


"""
    v = get_variate(support)

Compute a uniform variate within the `support`.
`support` is of the form `((x1,xspan),(y1,yspan),...)` where `xspan` is `x2 - x1`.
`x1` and `x2` are the bounds on `x`.

# Example
```jldoctest
julia> InterpolatedRejectionSampling.get_variate(((1,2),(10,1)))
(2.030427823129998, 10.800264489519098)
```
"""
function get_variate(support :: NTuple{N,NTuple{2,T}}) where {T,N}
    return ntuple( i -> (support[i][1] + rand()*support[i][2]), Val(N))
end

"""
univariate form
"""
function get_variate(support :: NTuple{2,T}) where {T}
    return support[1] + rand()*support[2]
end

"""
    v = get_variate(slice, support)

Compute a uniform variate within the sliced `support`.
`slice` is of the form `(x,y,:,...)` or `(:,:,z,...)` (or any permutation).
`x` (`y` and/or `z,...`) is a value that is not to be drawn but draw the remaining dimensions.
`support` is of the form `((x1,xspan),(y1,yspan),...)` where `xspan` is `x2 - x1`.
`x1` and `x2` are the bounds on `x`.

# Example
```jldoctest
julia> support = ((1,1),(2,3),(4,19))
((1, 1), (2, 3), (4, 19))

julia> InterpolatedRejectionSampling.get_variate(support)
(1.2865397858116667, 2.375903248257513, 17.251377235272265)

julia> slice = (1.1,missing,missing)
(1.1, missing, missing)

julia> InterpolatedRejectionSampling.get_variate(slice,support)
(1.1, 3.14019677035666, 6.252348183361507)
```
"""
function get_variate( slice   :: NTuple{N,Any}
                    , support :: NTuple{N,NTuple{2,Real}}) where N
    return ntuple( i -> ismissing(slice[i]) ? (support[i][1] + rand()*support[i][2]) : slice[i]
                 , Val(N)
                 )
end

"""
Perform a rejection sample

Given a `Interpolations.Extrapolation{T,N,ITPT,IT,ET}` construction, the support of each
dimension, and pmax
"""
function rsample( interp  :: Interpolations.Extrapolation{T,N,ITPT,IT,ET}
                , support :: NTuple{N,NTuple{2,Real}}
                , pmax    :: Real
                ) where {T,N,ITPT,IT,ET}
    # draw a sample
    samp = get_variate(support)
    # determine a rejection value
    prob = rand()*pmax

    # if rejection value is greater than the interpolated value of the function at sample
    while prob > interp(samp...)
        # redraw
        samp = get_variate(support)
        prob = rand()*pmax
    end
    return samp
end

"""
univariate form
"""
function rsample( interp  :: Interpolations.Extrapolation{T,1,ITPT,IT,ET}
                , support :: NTuple{2,Real}
                , pmax    :: Real
                ) where {T,ITPT,IT,ET}
    # draw a sample
    samp = get_variate(support)
    # determine a rejection value
    prob = rand()*pmax

    # if rejection value is greater than the interpolated value of the function at sample
    while prob > interp(samp)
        # redraw
        samp = get_variate(support)
        prob = rand()*pmax
    end
    return samp
end

"""
Perform a rejection sample for a slice

Given a `Interpolations.Extrapolation{T,N,ITPT,IT,ET}` construction, the support of each
dimension, and pmax
"""
function rsample( slice   :: NTuple{N,Any}
                , interp  :: Interpolations.Extrapolation{T,N,ITPT,IT,ET}
                , support :: NTuple{N,NTuple{2,Real}}
                , pmax :: Real
                ) where {T,N,ITPT,IT,ET}

    # draw a sample
    samp = get_variate(slice, support)
    # determine a rejection value
    prob = rand()*pmax

    # if rejection value is greater than the interpolated value of the function at sample
    while prob > interp(samp...)
        # redraw
        samp = get_variate(slice, support)
        prob = rand()*pmax
    end
    return samp
end

"""
    samples = irsample(knots, prob, n)

Draw `n` samples according to `prob` where `knots` are the dimensions along which `prob`
is constructed.

# Example
```jldoctest
julia> nx,ny = 10,5
julia> x_knots = range(-π/2,π/2,length=nx)
julia> y_knots = range(0,π,length=ny)
julia> prob = ones(nx,ny)
julia> for (i,p) in enumerate(cos.(x_knots))
       prob[i,:] .*= p
       end
julia> for (i,p) in enumerate(sin.(y_knots))
       prob[:,i] .*= p
       end
julia> irsample((x_knots, y_knots), prob, 5)
5-element Array{Tuple{Float64,Float64},1}:
 (1.3923419068184772, 1.654197559276566)
 (0.44884616227197993, 2.079250950153405)
 (0.007881035419223359, 1.4648986414076788)
 (-0.31237274389043135, 1.4515883747651241)
 (-1.2699084731307861, 0.660071653362384)
```
"""
function irsample( support :: NTuple{N,NTuple{2,Real}}
                 , interp  :: Interpolations.Extrapolation{T,N,ITPT,IT,ET}
                 , n       :: Integer
                 , pmax    :: Real
                 ) where {T,N,ITPT,IT,ET}
    return [rsample(interp,support,pmax) for i in 1:n]
end

"""
univariate form
"""
function irsample( support :: NTuple{2,Real}
                 , interp  :: Interpolations.Extrapolation{T,1,ITPT,IT,ET}
                 , n       :: Integer
                 , pmax    :: Real
                 ) where {T,ITPT,IT,ET}
    return [rsample(interp,support,pmax) for i in 1:n]
end

function irsample( knots :: NTuple{N,AbstractRange}
                 , prob  :: AbstractArray{T,N}
                 , n     :: Integer
                 , ::Linear
                 ; pmax  :: Real = PROB_PAD*maximum(prob)
                 ) where {T,N}
    return irsample(get_support(knots), LinearInterpolation(knots,prob),n,pmax)
end

function irsample( knots :: NTuple{N,AbstractRange}
                 , prob  :: AbstractArray{T,N}
                 , n     :: Integer
                 , ::CubicSpline
                 ; pmax  :: Real = PROB_PAD*maximum(prob)
                 ) where {T,N}
    return irsample(get_support(knots), CubicSplineInterpolation(knots,prob),n,pmax)
end

"""
univariate form
"""
function irsample( knots :: AbstractRange
                 , prob  :: AbstractVector
                 , n     :: Integer
                 , ::Linear
                 ; pmax  :: Real = PROB_PAD*maximum(prob)
                 )
    return irsample(get_support(knots), LinearInterpolation(knots,prob),n,pmax)
end

"""
univariate form
"""
function irsample( knots :: AbstractRange
                 , prob  :: AbstractVector
                 , n     :: Integer
                 , ::CubicSpline
                 ; pmax  :: Real = PROB_PAD*maximum(prob)
                 )
    return irsample(get_support(knots), CubicSplineInterpolation(knots,prob),n,pmax)
end

"""
default to CubicSpline when all knots are ranges
"""
function irsample( knots :: NTuple{N,AbstractRange}
                 , prob  :: AbstractArray{T,N}
                 , n     :: Integer
                 ; kwargs...
                 ) where {T,N}
    return irsample(knots,prob,n,CubicSpline();kwargs...)
end

"""
default to LinearInterpolation when knots are vectors
"""
function irsample( knots :: NTuple{N,AbstractVector}
                 , prob  :: AbstractArray{T,N}
                 , n     :: Integer
                 ; pmax  :: Real = PROB_PAD*maximum(prob)
                 ) where {T,N}
    return irsample(get_support(knots), LinearInterpolation(knots,prob),n,pmax)
end

"""
default to CubicSpline when all knots are ranges
univariate form
"""
function irsample( knots :: AbstractRange
                 , prob  :: AbstractVector
                 , n     :: Integer
                 ; kwargs...
                 )
    return irsample(knots,prob,n,CubicSpline();kwargs...)
end

"""
default to LinearInterpolation when knots are vectors
univariate form
"""
function irsample( knots :: AbstractVector
                 , prob  :: AbstractVector
                 , n     :: Integer
                 ; pmax  :: Real = PROB_PAD*maximum(prob)
                 )
    return irsample(get_support(knots), LinearInterpolation(knots,prob),n,pmax)
end

"""
    rejection_sampling!(slices, knots, prob)

`slices` is a `Vector` of `Tuples` of the form `(x,:,z,...)` where `:` is the dimension(s) to draw for the given values.
Samples are drawn according to `prob` where `knots` are used to construct the approximate `prob`
`size(prob) == length.(knots)`

# Example
```jldoctest
julia> nx,ny = 10,5
julia> x_knots = range(-π/2,π/2,length=nx)
julia> y_knots = range(0,π,length=ny)
julia> prob = ones(nx,ny)
julia> for (i,p) in enumerate(cos.(x_knots))
       prob[i,:] .*= p
       end
julia> for (i,p) in enumerate(sin.(y_knots))
       prob[:,i] .*= p
       end
julia> samples = [(0.1,:),(:,0.3)]
julia> rejection_sampling!(samples, (x_knots, y_knots), prob)
julia> display(samples)
2-element Array{Tuple{Float64,Float64},1}:
 (0.1, 2.0542954076141076)
 (-0.36452659968061485, 0.3)
```
"""
function irsample!( slices  :: Union{Vector,SubArray}
                  , support :: NTuple{N,NTuple{2,Real}}
                  , interp  :: Interpolations.Extrapolation{T,N,ITPT,IT,ET}
                  , pmax    :: Real
                  ) where {T,N,ITPT,IT,ET}
    for (i,slice) in enumerate(slices)
        slices[i] = rsample(slice,interp,support,pmax)
    end
end

function irsample!( slices :: Union{Vector,SubArray}
                  , knots  :: NTuple{N,AbstractRange}
                  , prob   :: AbstractArray{T,N}
                  , ::Linear
                  ; pmax   :: Real = PROB_PAD*maximum(prob)
                  ) where {T,N}
    irsample!(slices, get_support(knots), LinearInterpolation(knots,prob),pmax)
end

function irsample!( slices :: Union{Vector,SubArray}
                  , knots  :: NTuple{N,AbstractRange}
                  , prob   :: AbstractArray{T,N}
                  , ::CubicSpline
                  ; pmax   :: Real = PROB_PAD*maximum(prob)
                  ) where {T,N}
    irsample!(slices, get_support(knots), CubicSplineInterpolation(knots,prob),pmax)
end

"""
default to CubicSpline when all knots are ranges
"""
function irsample!( slices :: Union{Vector,SubArray}
                  , knots  :: NTuple{N,AbstractRange}
                  , prob  :: AbstractArray{T,N}
                  ; kwargs...
                  ) where {T,N}
    irsample!(slices,knots,prob,CubicSpline();kwargs...)
end

"""
default to LinearInterpolation when knots are not all ranges
"""
function irsample!( slices :: Union{Vector,SubArray}
                  , knots  :: NTuple{N,AbstractVector}
                  , prob   :: AbstractArray{T,N}
                  ; pmax  :: Real = PROB_PAD*maximum(prob)
                  ) where {T,N}
    irsample!(slices, get_support(knots), LinearInterpolation(knots,prob), pmax)
end

end
