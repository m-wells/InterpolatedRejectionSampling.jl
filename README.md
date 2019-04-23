# InterpolatedRejectionSampling.jl
[![Build Status](https://travis-ci.com/m-wells/InterpolatedRejectionSampling.jl.svg?token=qtRCxXQJn8B2HN1f6h3k&branch=master)](https://travis-ci.com/m-wells/InterpolatedRejectionSampling.jl)

## Draw samples from discrete multivariate distributions
For a given discrete (n-dimensional) grid of values (essentially the weights or probabilities) and the vectors that describe the span of the underlying space we can draw samples.
```julia
knots :: NTuple{N,V<:AbstractVector}
prob  :: AbstractArray{T<:Real,N}
```
The interpolation of the space is handled by  [`Interpolations.jl`](https://github.com/JuliaMath/Interpolations.jl)
## A simple example
First we need to setup a discrete distribution
```julia
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
```
Now we can call `irsample` (short form Interpolated Rejection Sample)
```
julia> irsample((x_knots, y_knots), prob, 5)
5-element Array{Tuple{Float64,Float64},1}:
 (1.3923419068184772, 1.654197559276566)
 (0.44884616227197993, 2.079250950153405)
 (0.007881035419223359, 1.4648986414076788)
 (-0.31237274389043135, 1.4515883747651241)
 (-1.2699084731307861, 0.660071653362384)
```
`irsample` returns an array of `Tuples` because this in a n-dimensional sampling we are drawing points, hence `Tuples`

## More complicated
Go look at the `test/runtests.jl` file
