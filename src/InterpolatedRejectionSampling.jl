#=
    InterpolatedRejectionSampling
    Copyright © 2019 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#

module InterpolatedRejectionSampling

using Interpolations

export rejection_sampling

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

julia> slice = (1.1,:,:)
(1.1, Colon(), Colon())

julia> InterpolatedRejectionSampling.get_variate(slice,support)
(1.1, 3.14019677035666, 6.252348183361507)
```
"""
function get_variate( slice   :: NTuple{N,Union{Real,Colon}}
                    , support :: NTuple{N,NTuple{2,Real}}) where N
    return ntuple( i -> isa(slice[i],Colon) ? (support[i][1] + rand()*support[i][2]) : slice[i]
                 , Val(N)
                 )
end

function rsample( sitp    :: Interpolations.Extrapolation{T,N,ITPT,IT,ET}
                , support :: NTuple{N,NTuple{2,Real}}
                , pmax    :: Real
                ) where {T,N,ITPT,IT,ET}
    # draw a sample
    samp = get_variate(support)
    # determine a rejection value
    prob = rand()*pmax

    # if rejection value is greater than the interpolated value of the function at sample
    while prob > sitp(samp...)
        # redraw
        samp = get_variate(support)
        prob = rand()*pmax
    end
    return samp
end

function rsample( sitp    :: Interpolations.Extrapolation{T,N,ITPT,IT,ET}
                , slice   :: NTuple{N,Union{Real,Colon}}
                , support :: NTuple{N,NTuple{2,Real}}
                , pmax :: Real
                ) where {T,N,ITPT,IT,ET}

    # draw a sample
    samp = get_variate(slice, support)
    # determine a rejection value
    prob = rand()*pmax

    # if rejection value is greater than the interpolated value of the function at sample
    while prob > sitp(samp...)
        # redraw
        samp = get_variate(slice, support)
        prob = rand()*pmax
    end
    return samp
end

"""
    samples = rejection_sampling(knots, prob, n)

Draw `n` samples according to `prob` where `knots` where used to construct the approximate `prob`

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
julia> rejection_sampling((x_knots, y_knots), prob, 5)
5-element Array{Tuple{Float64,Float64},1}:
 (1.3923419068184772, 1.654197559276566)
 (0.44884616227197993, 2.079250950153405)
 (0.007881035419223359, 1.4648986414076788)
 (-0.31237274389043135, 1.4515883747651241)
 (-1.2699084731307861, 0.660071653362384)
```
"""
function rejection_sampling( knots :: NTuple{N,AbstractRange}
                           , prob  :: AbstractArray{T,N}
                           , n     :: Integer
                           ) where {T,N}
    # interpolate the prob with cubic spline with linear edges with knots on cell edges
    cbs_interp = CubicSplineInterpolation(knots, prob)
    support = ntuple( i -> ( knots[i][1]
                           , knots[i][end] - knots[i][1]
                           )
                    , Val(N)
                    )

    samp = Vector{NTuple{N,T}}(undef,n)
    pmax = 1.05*maximum(prob)
    for i in eachindex(samp)
        samp[i] = rsample(cbs_interp,support,pmax)
    end

    return samp
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
function rejection_sampling!( slices :: AbstractVector{NTuple{N,Any}}
                            , knots  :: NTuple{N,AbstractRange}
                            , prob   :: AbstractArray{T,N}
                           ) where {N,T}
    # interpolate the prob with cubic spline with linear edges with knots on cell edges
    cbs_interp = CubicSplineInterpolation(knots, prob)
    support = ntuple( i -> ( knots[i][1]
                           , knots[i][end] - knots[i][1]
                           )
                    , Val(N)
                    )

    pmax = 1.05*maximum(prob)
    for (i,slice) in enumerate(slices)
        slices[i] = rsample(cbs_interp,slice,support,pmax)
    end

    return nothing
end
####################################################################################################
#
####################################################################################################
#---------------------------------------------------------------------------------------------------
end
