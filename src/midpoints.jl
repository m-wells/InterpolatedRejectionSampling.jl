# Very efficient midpoint calculators

@inline function midpoints(x::AbstractVector{T}) where T<:Real
    retval = Vector{Float64}(undef, length(x)-1)
    @fastmath @inbounds @simd for i in eachindex(retval)
        retval[i] = (x[i  ] +
                     x[i+1])/2
    end
    return retval
end

@inline function midpoints(x::AbstractArray{T,2}) where T<:Real
    retval = Matrix{Float64}(undef,size(x).-1)
    @fastmath @inbounds for j in 1:last(size(retval))
    @simd for i in 1:first(size(retval))
        retval[i,j] = (x[i  ,j  ] +
                       x[i+1,j  ] +
                       x[i  ,j+1] +
                       x[i+1,j+1])/4
    end; end
    return retval
end

@inline function midpoints(x::AbstractArray{T,3}) where T<:Real
    retval = Array{Float64,3}(undef,size(x).-1)
    @fastmath @inbounds for k in 1:last(size(retval))
    for j in 1:size(retval)[2]
    @simd for i in 1:first(size(retval))
        retval[i,j,k] = (x[i  ,j  ,k  ] +
                         x[i+1,j  ,k  ] +
                         x[i  ,j+1,k  ] +
                         x[i+1,j+1,k  ] +
                         x[i  ,j  ,k+1] +
                         x[i+1,j  ,k+1] +
                         x[i  ,j+1,k+1] +
                         x[i+1,j+1,k+1])/8
    end; end; end
    return retval
end

@inline function midpoints(x::AbstractArray{T,4}) where T<:Real
    retval = Array{Float64,4}(undef,size(x).-1)
    @fastmath @inbounds for l in 1:last(size(retval))
    for k in 1:size(retval)[3]
    for j in 1:size(retval)[2]
    @simd for i in 1:first(size(retval))
        retval[i,j,k,l] = (x[i  ,j  ,k  ,l  ] +
                           x[i+1,j  ,k  ,l  ] +
                           x[i  ,j+1,k  ,l  ] +
                           x[i+1,j+1,k  ,l  ] +
                           x[i  ,j  ,k+1,l  ] +
                           x[i+1,j  ,k+1,l  ] +
                           x[i  ,j+1,k+1,l  ] +
                           x[i+1,j+1,k+1,l  ] +
                           x[i  ,j  ,k  ,l+1] +
                           x[i+1,j  ,k  ,l+1] +
                           x[i  ,j+1,k  ,l+1] +
                           x[i+1,j+1,k  ,l+1] +
                           x[i  ,j  ,k+1,l+1] +
                           x[i+1,j  ,k+1,l+1] +
                           x[i  ,j+1,k+1,l+1] +
                           x[i+1,j+1,k+1,l+1])/16
    end; end; end; end
    return retval
end

@inline function midpoints(x::AbstractArray{T,5}) where T<:Real
    retval = Array{Float64,5}(undef,size(x).-1)
    @fastmath @inbounds for m in 1:last(size(retval))
    for l in 1:size(retval)[4]
    for k in 1:size(retval)[3]
    for j in 1:size(retval)[2]
    @simd for i in 1:first(size(retval))
        retval[i,j,k,l,m] = (x[i  ,j  ,k  ,l  ,m  ] +
                             x[i+1,j  ,k  ,l  ,m  ] +
                             x[i  ,j+1,k  ,l  ,m  ] +
                             x[i+1,j+1,k  ,l  ,m  ] +
                             x[i  ,j  ,k+1,l  ,m  ] +
                             x[i+1,j  ,k+1,l  ,m  ] +
                             x[i  ,j+1,k+1,l  ,m  ] +
                             x[i+1,j+1,k+1,l  ,m  ] +
                             x[i  ,j  ,k  ,l+1,m  ] +
                             x[i+1,j  ,k  ,l+1,m  ] +
                             x[i  ,j+1,k  ,l+1,m  ] +
                             x[i+1,j+1,k  ,l+1,m  ] +
                             x[i  ,j  ,k+1,l+1,m  ] +
                             x[i+1,j  ,k+1,l+1,m  ] +
                             x[i  ,j+1,k+1,l+1,m  ] +
                             x[i+1,j+1,k+1,l+1,m  ] +
                             x[i  ,j  ,k  ,l  ,m+1] +
                             x[i+1,j  ,k  ,l  ,m+1] +
                             x[i  ,j+1,k  ,l  ,m+1] +
                             x[i+1,j+1,k  ,l  ,m+1] +
                             x[i  ,j  ,k+1,l  ,m+1] +
                             x[i+1,j  ,k+1,l  ,m+1] +
                             x[i  ,j+1,k+1,l  ,m+1] +
                             x[i+1,j+1,k+1,l  ,m+1] +
                             x[i  ,j  ,k  ,l+1,m+1] +
                             x[i+1,j  ,k  ,l+1,m+1] +
                             x[i  ,j+1,k  ,l+1,m+1] +
                             x[i+1,j+1,k  ,l+1,m+1] +
                             x[i  ,j  ,k+1,l+1,m+1] +
                             x[i+1,j  ,k+1,l+1,m+1] +
                             x[i  ,j+1,k+1,l+1,m+1] +
                             x[i+1,j+1,k+1,l+1,m+1])/32
    end; end; end; end; end
    return retval
end

"""
the inner most loop needs to be broken into two additions
if not then it overflows the cache and allocations go crazy

midpoints(rand(2,2,2,2,2,2)) allocates 6 times with 384 bytes like this
when trying to shove all 64 onto a single line it had 2.15k allocations with 49.141 KiB
"""
@inline function midpoints(x::AbstractArray{T,6}) where T<:Real
    retval = Array{Float64,6}(undef,size(x).-1)
    @fastmath @inbounds for n in 1:last(size(retval))
    for m in 1:size(retval)[5]
    for l in 1:size(retval)[4]
    for k in 1:size(retval)[3]
    for j in 1:size(retval)[2]
    @simd for i in 1:first(size(retval))
        retval[i,j,k,l,m,n] = (x[i  ,j  ,k  ,l  ,m  ,n  ] +
                               x[i+1,j  ,k  ,l  ,m  ,n  ] +
                               x[i  ,j+1,k  ,l  ,m  ,n  ] +
                               x[i+1,j+1,k  ,l  ,m  ,n  ] +
                               x[i  ,j  ,k+1,l  ,m  ,n  ] +
                               x[i+1,j  ,k+1,l  ,m  ,n  ] +
                               x[i  ,j+1,k+1,l  ,m  ,n  ] +
                               x[i+1,j+1,k+1,l  ,m  ,n  ] +
                               x[i  ,j  ,k  ,l+1,m  ,n  ] +
                               x[i+1,j  ,k  ,l+1,m  ,n  ] +
                               x[i  ,j+1,k  ,l+1,m  ,n  ] +
                               x[i+1,j+1,k  ,l+1,m  ,n  ] +
                               x[i  ,j  ,k+1,l+1,m  ,n  ] +
                               x[i+1,j  ,k+1,l+1,m  ,n  ] +
                               x[i  ,j+1,k+1,l+1,m  ,n  ] +
                               x[i+1,j+1,k+1,l+1,m  ,n  ] +
                               x[i  ,j  ,k  ,l  ,m+1,n  ] +
                               x[i+1,j  ,k  ,l  ,m+1,n  ] +
                               x[i  ,j+1,k  ,l  ,m+1,n  ] +
                               x[i+1,j+1,k  ,l  ,m+1,n  ] +
                               x[i  ,j  ,k+1,l  ,m+1,n  ] +
                               x[i+1,j  ,k+1,l  ,m+1,n  ] +
                               x[i  ,j+1,k+1,l  ,m+1,n  ] +
                               x[i+1,j+1,k+1,l  ,m+1,n  ] +
                               x[i  ,j  ,k  ,l+1,m+1,n  ] +
                               x[i+1,j  ,k  ,l+1,m+1,n  ] +
                               x[i  ,j+1,k  ,l+1,m+1,n  ] +
                               x[i+1,j+1,k  ,l+1,m+1,n  ] +
                               x[i  ,j  ,k+1,l+1,m+1,n  ] +
                               x[i+1,j  ,k+1,l+1,m+1,n  ] +
                               x[i  ,j+1,k+1,l+1,m+1,n  ] +
                               x[i+1,j+1,k+1,l+1,m+1,n  ])/32
        retval[i,j,k,l,m,n] += (x[i  ,j  ,k  ,l  ,m  ,n+1] +
                                x[i+1,j  ,k  ,l  ,m  ,n+1] +
                                x[i  ,j+1,k  ,l  ,m  ,n+1] +
                                x[i+1,j+1,k  ,l  ,m  ,n+1] +
                                x[i  ,j  ,k+1,l  ,m  ,n+1] +
                                x[i+1,j  ,k+1,l  ,m  ,n+1] +
                                x[i  ,j+1,k+1,l  ,m  ,n+1] +
                                x[i+1,j+1,k+1,l  ,m  ,n+1] +
                                x[i  ,j  ,k  ,l+1,m  ,n+1] +
                                x[i+1,j  ,k  ,l+1,m  ,n+1] +
                                x[i  ,j+1,k  ,l+1,m  ,n+1] +
                                x[i+1,j+1,k  ,l+1,m  ,n+1] +
                                x[i  ,j  ,k+1,l+1,m  ,n+1] +
                                x[i+1,j  ,k+1,l+1,m  ,n+1] +
                                x[i  ,j+1,k+1,l+1,m  ,n+1] +
                                x[i+1,j+1,k+1,l+1,m  ,n+1] +
                                x[i  ,j  ,k  ,l  ,m+1,n+1] +
                                x[i+1,j  ,k  ,l  ,m+1,n+1] +
                                x[i  ,j+1,k  ,l  ,m+1,n+1] +
                                x[i+1,j+1,k  ,l  ,m+1,n+1] +
                                x[i  ,j  ,k+1,l  ,m+1,n+1] +
                                x[i+1,j  ,k+1,l  ,m+1,n+1] +
                                x[i  ,j+1,k+1,l  ,m+1,n+1] +
                                x[i+1,j+1,k+1,l  ,m+1,n+1] +
                                x[i  ,j  ,k  ,l+1,m+1,n+1] +
                                x[i+1,j  ,k  ,l+1,m+1,n+1] +
                                x[i  ,j+1,k  ,l+1,m+1,n+1] +
                                x[i+1,j+1,k  ,l+1,m+1,n+1] +
                                x[i  ,j  ,k+1,l+1,m+1,n+1] +
                                x[i+1,j  ,k+1,l+1,m+1,n+1] +
                                x[i  ,j+1,k+1,l+1,m+1,n+1] +
                                x[i+1,j+1,k+1,l+1,m+1,n+1])/32
    end; end; end; end; end; end
    return retval
end
