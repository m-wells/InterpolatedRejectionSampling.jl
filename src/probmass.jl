include("midpoints.jl")

@inline function probmass(x::AbstractVector{T}, y::AbstractVector{T}
                         ) where T<:Real
    retval = midpoints(y)
    Δx = diff(x)
    @fastmath @inbounds @simd for i in eachindex(Δx)
        retval[i] *= Δx[i]
    end
    return retval./sum(retval)
end

@inline function probmass(x::NTuple{2,AbstractVector{T}}, y::AbstractMatrix{T}
                         ) where T<:Real
    retval = midpoints(y)
    Δx1,Δx2 = diff.(x)
    @fastmath @inbounds for j in eachindex(Δx2)
    @simd for i in eachindex(Δx1)
        retval[i,j] *= Δx1[i]*Δx2[j]
    end; end
    return retval./sum(retval)
end

@inline function probmass(x::NTuple{3,AbstractVector{T}}, y::AbstractArray{T,3}
                         ) where T<:Real
    retval = midpoints(y)
    Δx1,Δx2,Δx3 = diff.(x)
    @fastmath @inbounds for k in eachindex(Δx3)
    for j in eachindex(Δx2)
    @simd for i in eachindex(Δx1)
        retval[i,j,k] *= Δx1[i]*Δx2[j]*Δx3[k]
    end; end; end
    return retval./sum(retval)
end

@inline function probmass(x::NTuple{4,AbstractVector{T}}, y::AbstractArray{T,4}
                         ) where T<:Real
    retval = midpoints(y)
    Δx1,Δx2,Δx3,Δx4 = diff.(x)
    @fastmath @inbounds for l in eachindex(Δx4)
    for k in eachindex(Δx3)
    for j in eachindex(Δx2)
    @simd for i in eachindex(Δx1)
        retval[i,j,k,l] *= Δx1[i]*Δx2[j]*Δx3[k]*Δx4[l]
    end; end; end; end
    return retval./sum(retval)
end

@inline function probmass(x::NTuple{5,AbstractVector{T}}, y::AbstractArray{T,5}
                         ) where T<:Real
    retval = midpoints(y)
    Δx1,Δx2,Δx3,Δx4,Δx5 = diff.(x)
    @fastmath @inbounds for m in eachindex(Δx5)
    for l in eachindex(Δx4)
    for k in eachindex(Δx3)
    for j in eachindex(Δx2)
    @simd for i in eachindex(Δx1)
        retval[i,j,k,l,m] *= Δx1[i]*Δx2[j]*Δx3[k]*Δx4[l]*Δx5[m]
    end; end; end; end; end
    return retval./sum(retval)
end

@inline function probmass(x::NTuple{6,AbstractVector{T}}, y::AbstractArray{T,6}
                         ) where T<:Real
    retval = midpoints(y)
    Δx1,Δx2,Δx3,Δx4,Δx5,Δx6 = diff.(x)
    @fastmath @inbounds for n in eachindex(Δx6)
    for m in eachindex(Δx5)
    for l in eachindex(Δx4)
    for k in eachindex(Δx3)
    for j in eachindex(Δx2)
    @simd for i in eachindex(Δx1)
        retval[i,j,k,l,m,n] *= Δx1[i]*Δx2[j]*Δx3[k]*Δx4[l]*Δx5[m]*Δx6[n]
    end; end; end; end; end; end
    return retval./sum(retval)
end
