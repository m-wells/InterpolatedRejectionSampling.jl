#using Interpolations: GridIndex, tweight, tcoef, check_gridded
#
#mutable struct MutableExtrapolation{T,N,ITPT,IT,ET} <: AbstractExtrapolation{T,N,ITPT,IT}
#    itp::ITPT
#    et::ET
#end
#
#mutable struct MutableGriddedInterpolation{T,N,TCoefs,IT<:DimSpec{Gridded},K<:Tuple{Vararg{AbstractVector}}} <: AbstractInterpolation{T,N,IT}
#    knots::K
#    coefs::Array{TCoefs,N}
#    it::IT
#end
#
#lbounds(itp::MutableGriddedInterpolation) = first.(itp.knots)
#ubounds(itp::MutableGriddedInterpolation) = last.(itp.knots)
#
#
#function MutableGriddedInterpolation(::Type{TWeights}, knots::NTuple{N,GridIndex}, A::AbstractArray{TCoefs,N}, it::IT) where {N,TCoefs,TWeights<:Real,IT<:DimSpec{Gridded},pad}
#    isconcretetype(IT) || error("The b-spline type must be a leaf type (was $IT)")
#    isconcretetype(TCoefs) || @warn("For performance reasons, consider using an array of a concrete type (eltype(A) == $(eltype(A)))")
#
#    check_gridded(it, knots, axes(A))
#    c = zero(TWeights)
#    if isempty(A)
#        T = Base.promote_op(*, typeof(c), eltype(A))
#    else
#        T = typeof(c * first(A))
#    end
#    MutableGriddedInterpolation{T,N,TCoefs,IT,typeof(knots)}(knots, A, it)
#end
#
#function mutable_interpolate(::Type{TWeights}, ::Type{TCoefs}, knots::NTuple{N,GridIndex}, A::AbstractArray{Tel,N}, it::IT) where {TWeights,TCoefs,Tel,N,IT<:DimSpec{Gridded}}
#    MutableGriddedInterpolation(TWeights, knots, A, it)
#end
#
#function mutable_interpolate(knots::NTuple{N,GridIndex}, A::AbstractArray{Tel,N}, it::IT) where {Tel,N,IT<:DimSpec{Gridded}}
#    mutable_interpolate(tweight(A), tcoef(A), knots, A, it)
#end
#
#MutableLinearInterpolation(ranges::NTuple{N,AbstractVector},
#                           vs::AbstractArray{T,N};
#                           extrapolation_bc = Throw()
#                          ) where {N,T} = extrapolate(mutable_interpolate(ranges, vs, Gridded(Linear())), extrapolation_bc)
