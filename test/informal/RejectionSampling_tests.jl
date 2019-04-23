#=
    RejectionSampling_tests
    Copyright © 2019 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#
using StatsBase
using PyPlot
import LinearAlgebra.normalize!
using Distributions
using PyCall
axes_grid1 = pyimport("mpl_toolkits.axes_grid1")



include("../src/InterpolatedRejectionSampling.jl")
####################################################################################################
function test_dist1(x,y)
    prob_matrix = Matrix{Float64}(undef,length(x),length(y))
    probx = [pdf(Beta(2,5), t) for t in x]
    proby = [pdf(Beta(1,3), t) for t in y]

    for (j,py) in enumerate(proby)
        for (i,px) in enumerate(probx)
            prob_matrix[i,j] = px*py
        end
    end
    return prob_matrix
end
x1,x2 = 0,1
y1,y2 = 0,1

x = range(x1,x2,length=10)
y = range(y1,y2,length=5)
prob_matrix = test_dist1(x,y)
display(prob_matrix)

samps = InterpolatedRejectionSampling.rejection_sampling((x,y), prob_matrix, 1000000)
xs = getindex.(samps,1)
ys = getindex.(samps,2)

hist = fit( Histogram
          , (xs,ys)
          , ( range(x1,x2,length=100)
            , range(y1,y2,length=100)
            )
          ; closed = :left
          )

axes = axes_grid1[:ImageGrid](plt[:figure](), 111,
                 nrows_ncols=(1,2),
                 axes_pad=0.3,
                 share_all=true,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="2%",
                 cbar_pad=0.15,
                 label_mode = "L"
                 )

img = axes[1][:imshow]( hist.weights
                , origin = "lower"
                , extent = (x1,x2,y1,y2)
                )

x = range(x1,x2,length=1000)
y = range(y1,y2,length=1000)
prob_matrix = test_dist1(x,y)
#X = [i for i=eachindex(x),j=eachindex(y)]
#Y = [j for i=eachindex(x),j=eachindex(y)]

mask = @. !isinf(prob_matrix) & !isnan(prob_matrix)
prob_matrix ./= sum(prob_matrix[mask])
axes[2][:contourf]( prob_matrix
                  , 20
                  , origin = "lower"
                  , extent = (x1,x2,y1,y2)
                  )

axes[end][:cax][:colorbar](img)

####################################################################################################
#
####################################################################################################
#axes[2][:contour](X,Y,prob_matrix)

#@code_warntype RejectionSampling.rejection_sampling(5, pdf, x, y, z)
#@code_warntype RejectionSampling.rejection_sampling(5, pdf, x, y)

#function array3(fillval, N)
#    fill(fillval, ntuple(d->3, N))
#end
#
#array3(5.0, 2)
#@code_warntype array3(5.0, 2)
