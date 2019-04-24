#=
    run_tests
    Copyright © 2019 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the MIT license.
=#

using InterpolatedRejectionSampling
using Test

import Random.seed!
seed!(1234)

function scale_range( x :: AbstractVector
                    , a :: Real
                    , b :: Real
                    )   :: Vector{Float64}
    return (b-a)*(x .- minimum(x))/(maximum(x) - minimum(x)) .+ a
end

function mean( x :: AbstractVector )
    return sum(x)/length(x)
end

@testset "rejection sample" begin
    x1,x2 = 0,π
    y1,y2 = -π/2,π/2

    X = range(x1,x2,length=20)
    Y = scale_range(sort(rand(10)),y1,y2)
    knots = (X,Y)

    fx(x) = sin(x)
    fy(y) = cos(y)

    A = Array{Float64,2}(undef,length(X),length(Y))
    for (i,x) in enumerate(X)
        for (j,y) in enumerate(Y)
            A[i,j] = fx(x)+fy(y)
        end
    end
    s = irsample(knots,A,100000)
    sx = [xy[1] for xy in s]
    sy = [xy[2] for xy in s]
    @test isapprox(mean(sx),π/2,atol=1e-2)
    @test isapprox(mean(sy),0.0,atol=1e-2)
end

@testset "sliced rejection sample" begin
    x1,x2 = 0,π
    y1,y2 = -π/2,π/2
    z1,z2 = -π/4,π/4

    X = collect(range(x1,x2,length=20))
    Y = collect(range(y1,y2,length=10))
    Z = scale_range(sort(rand(10)),z1,z2)
    knots = (X,Y,Z)

    fx(x) = sin(x)
    fy(y) = cos(y)
    fz(z) = tan(z)

    A = Array{Float64,3}(undef,length(X),length(Y),length(Z))
    for (i,x) in enumerate(X)
        for (j,y) in enumerate(Y)
            for (k,z) in enumerate(Z)
                A[i,j,k] = fx(x)+fy(y)+fz(z)
            end
        end
    end

    a = 0.1
    b = π/2
    c = -1.0
    d = 0.0

    sx = [a,b,missing,missing,a]
    sy = [missing,missing,c,missing,b]
    sz = [missing,d,d,missing,d]
    szip = zip(sx,sy,sz)

    # Zip, Vector, SubArray
    for s in [collect(szip),selectdim(collect(szip),1,:)]
        irsample!(s,knots,A)

        @test s[1][1] == a
        @test y1 ≤ s[1][2] ≤ y2
        @test z1 ≤ s[1][3] ≤ z2

        @test s[2][1] == b
        @test y1 ≤ s[2][2] ≤ y2
        @test s[2][3] == d

        @test x1 ≤ s[3][1] ≤ x2
        @test s[3][2] == c
        @test s[3][3] == d

        @test x1 ≤ s[4][1] ≤ x2
        @test y1 ≤ s[4][2] ≤ y2
        @test z1 ≤ s[4][3] ≤ z2

        @test s[5][1] == a
        @test s[5][2] == b
        @test s[5][3] == d
    end
end

@testset "univariate sample" begin
    x1,x2 = 0,π
    x = range(x1,x2,length=20)
    y = sin.(x)
    s = irsample(x,y,100000)
    @test isapprox(mean(s),π/2,atol=1e-2)
end
