#=
    slice_test
    Copyright © 2019 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the MIT license.
=#

include("../src/InterpolatedRejectionSampling.jl")

if false
    nx,ny = 10,5
    x_knots = range(-π/2,π/2,length=nx)
    y_knots = range(0,π,length=ny)
    
    prob = ones(nx,ny)
    for (i,p) in enumerate(cos.(x_knots))
        prob[i,:] .*= p
    end
    for (i,p) in enumerate(sin.(y_knots))
        prob[:,i] .*= p
    end
    
    InterpolatedRejectionSampling.rejection_sampling(5, prob, x_knots, y_knots)
end

#function somefunction1()
#    nx,ny = 10,5
#    x_knots = range(-π/2,π/2,length=nx)
#    y_knots = range(0,π,length=ny)
#    
#    prob = ones(nx,ny)
#    for (i,p) in enumerate(cos.(x_knots))
#        prob[i,:] .*= p
#    end
#    for (i,p) in enumerate(sin.(y_knots))
#        prob[:,i] .*= p
#    end
#
#    samples = InterpolatedRejectionSampling.rejection_sampling( (x_knots, y_knots)
#                                                          , prob
#                                                          , 5
#                                                          )
#    display(samples)
#end
#@code_warntype somefunction1()


#function somefunction2()
#    nx,ny = 10,5
#    x_knots = range(-π/2,π/2,length=nx)
#    y_knots = range(0,π,length=ny)
#    
#    prob = ones(nx,ny)
#    for (i,p) in enumerate(cos.(x_knots))
#        prob[i,:] .*= p
#    end
#    for (i,p) in enumerate(sin.(y_knots))
#        prob[:,i] .*= p
#    end
#
#    samples = [(0.1,:),(:,0.3)]
#    InterpolatedRejectionSampling.rejection_sampling!( samples
#                                                     , (x_knots, y_knots)
#                                                     , prob
#                                                     )
#    display(samples)
#end
#@code_warntype somefunction2()
