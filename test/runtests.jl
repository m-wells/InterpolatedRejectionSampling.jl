#=
    run_tests
    Copyright © 2019 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the MIT license.
=#

using InterpolatedRejectionSampling
using Test
using Statistics

import Random.seed!
seed!(1234)

function scale_range( x :: AbstractVector
                    , a :: Real
                    , b :: Real
                    )   :: Vector{Float64}
    return (b-a)*(x .- minimum(x))/(maximum(x) - minimum(x)) .+ a
end

@testset "rejection sample" begin
    fx(x) = sin(x)
    fy(y) = cos(y)

    X = range(0,π,length=20)
    Y = scale_range(sort(rand(10)),-π/2,π/2)
    A = Array{Float64,2}(undef,length(X),length(Y))
    for (i,x) in enumerate(X)
        for (j,y) in enumerate(Y)
            A[i,j] = fx(x)+fy(y)
        end
    end
    knots = (X,Y)
    s = irsample(knots,A,100000)
    sx = [xy[1] for xy in s]
    sy = [xy[2] for xy in s]
    @test isapprox(mean(sx),π/2,atol=1e-2)
    @test isapprox(mean(sy),0.0,atol=1e-2)
end

@testset "sliced rejection sample" begin
    fx(x) = sin(x)
    fy(y) = cos(y)
    fz(z) = tan(z)

    X = range(0,π,length=20)
    Y = range(-π/2,π/2,length=10)
    Z = scale_range(sort(rand(10)),-π/2,π/2)
    A = Array{Float64,2}(undef,length(X),length(Y),length(Z))
    for (i,x) in enumerate(X)
        for (j,y) in enumerate(Y)
            for (k,z) in enumerate(Z)
                A[i,j,k] = fx(x)+fy(y)+fz(z)
            end
        end
    end
    knots = (X,Y,Z)
    s = irsample(knots,A,100000)
    sx = [xy[1] for xy in s]
    sy = [xy[2] for xy in s]
    @test isapprox(mean(sx),π/2,atol=1e-2)
    @test isapprox(mean(sy),0.0,atol=1e-2)
end
