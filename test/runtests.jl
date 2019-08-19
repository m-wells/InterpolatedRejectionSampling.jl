using InterpolatedRejectionSampling
using Test

import Random.seed!
seed!(1234)

@testset "irsample" begin
    X = range(0.0, float(π), length=5)
    Y = range(0.0, π/4, length=4)
    n = 100

    x = irsample(X, sin.(X), n)
    @test isa(x, Vector{Float64})
    @test length(x) == n

    knots = (X,Y)
    prob = [sin(x)+tan(y) for x in X, y in Y]
    
    xy = irsample(knots, prob, n)
    @test isa(xy, Matrix{Float64})
    @test size(xy) == (2,n)

    xy = Matrix{Union{Float64,Missing}}(missing,2,n)
    irsample!(xy, knots, prob)
    @test isa(xy, Matrix{Union{Missing,Float64}})
    @test size(xy) == (2,n)
    @test iszero(count(ismissing.(xy)))

    xy = Matrix{Union{Float64,Missing}}(missing,2,n)
    xy[1,:] .= π.*rand(n)
    irsample!(xy, knots, prob)
    @test isa(xy, Matrix{Union{Missing,Float64}})
    @test size(xy) == (2,n)
    @test iszero(count(ismissing.(xy)))

    xy = Matrix{Union{Float64,Missing}}(missing,2,n)
    xy[1,1:2:n] .= π.*rand(length(1:2:n))
    irsample!(xy, knots, prob)
    @test isa(xy, Matrix{Union{Missing,Float64}})
    @test size(xy) == (2,n)
    @test iszero(count(ismissing.(xy)))
end 
