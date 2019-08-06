#=
    run_tests
    Copyright © 2019 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the MIT license.
=#

using InterpolatedRejectionSampling
using Test

import Random.seed!
seed!(1234)

@testset "irsample" begin

    X = range(0.0;stop=float(π),length=5)
    Y = range(0.0;stop=π/4,length=4)
    knots = (X,Y)
    prob = [sin(x)+tan(y) for x in X, y in Y]
    
    n = 100
    xy = irsample(knots, prob, n)
    @test isa(xy, Matrix{Float64})
    @test size(xy) == (2,n)

    xy = Matrix{Union{Float64,Missing}}(missing,2,n)
    irsample!(xy, knots, prob)
    @test isa(xy, Matrix{Union{Missing,Float64}})
    @test size(xy) == (2,n)
    @test iszero(count(ismissing.(xy)))

    x = irsample(X, sin.(X), n)
    @test isa(x, Vector{Float64})
    @test length(x) == n
end 
