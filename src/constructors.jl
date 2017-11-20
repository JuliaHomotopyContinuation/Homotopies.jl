# using Homotopy
# using DynamicPolynomials: @polyvar
# import MultivariatePolynomials
# const MP = MultivariatePolynomials
# import FixedPolynomials
# const FP = FixedPolynomials
#
# @polyvar x y z
#
# f = FP.Polynomial(x^2+3y+z)
# g = FP.Polynomial(x^3+z+2.0y^2)
#
#

function construct(::Type{AH},st::Tuple{Vector{FP.Polynomial{T}}, Vector{FP.Polynomial{T}}}
    ) where {T<:Number, AH<:AbstractPolynomialHomotopy}
    start, target = st
    AH{T}(start, target)
end

function construct(::Type{AH},st::Tuple{Vector{FP.Polynomial{U}}, Vector{FP.Polynomial{V}}}
    ) where {T<:Number, AH<:AbstractPolynomialHomotopy{T}, U<:Number, V<:Number}
    start, target = st
    AH(convert(Vector{FP.Polynomial{T}}, start), convert(Vector{FP.Polynomial{T}}, target))
end

function construct(::Type{AH},st::Tuple{Vector{FP.Polynomial{U}}, Vector{FP.Polynomial{V}}}
    ) where {AH<:AbstractPolynomialHomotopy, U<:Number, V<:Number}
    start, target = st
    P = promote_type(U, V, Float64)
    AH{P}(convert(Vector{FP.Polynomial{P}}, start),
        convert(Vector{FP.Polynomial{P}}, target))
end

function construct(::Type{AH},st::Tuple{Vector{<:MP.AbstractPolynomial{U}}, Vector{<:MP.AbstractPolynomial{V}}}
    ) where {T<:Number, AH<:AbstractPolynomialHomotopy{T}, U<:Number, V<:Number}
    start, target = st
    AH(convert(Vector{FP.Polynomial{T}}, start), convert(Vector{FP.Polynomial{T}}, target))
end

function construct(::Type{AH},st::Tuple{Vector{<:MP.AbstractPolynomial{U}}, Vector{<:MP.AbstractPolynomial{V}}}
    ) where {AH<:AbstractPolynomialHomotopy, U<:Number, V<:Number}
    start, target = st
    P = promote_type(U, V, Float64)
    AH(convert(Vector{FP.Polynomial{P}}, start),
        convert(Vector{FP.Polynomial{P}}, target))
end

function construct(::Type{AH},st::Tuple{<:MP.AbstractPolynomial{U}, <:MP.AbstractPolynomial{V}}
    ) where {T<:Number, AH<:AbstractPolynomialHomotopy{T}, U<:Number, V<:Number}
    start, target = st
    s, t = convert(Vector{FP.Polynomial{T}}, [start, target])
    AH([s], [t])
end

function construct(::Type{AH},st::Tuple{<:MP.AbstractPolynomial{U}, <:MP.AbstractPolynomial{V}}
    ) where {AH<:AbstractPolynomialHomotopy, U<:Number, V<:Number}
    start, target = st
    P = promote_type(U, V, Float64)
    start, target = st
    s, t = convert(Vector{FP.Polynomial{P}}, [start, target])
    AH{P}([s], [t])
end
#
#
# H = convert(StraightLineHomotopy{Complex128}, ([f], [g]))
# typeof(H)
#
# H = convert(StraightLineHomotopy, ([f], [g]))
#
# H = convert(StraightLineHomotopy, (2+3x+y, 2x+3.0y))
# H = convert(GeodesicOnTheSphere{Complex128}, (2+3x+y, 2x+3.0y))
# typeof(H)
#
#
#
#
# function GeodesicOnTheSphere{T}(
#     start::Vector{FP.Polynomial{U}},
#     target::Vector{FP.Polynomial{V}}
#     ) where {T<:Number, U<:Number, V<:Number}
#     GeodesicOnTheSphere{T}(
#         convert(Vector{FP.Polynomial{T}}, start),
#         convert(Vector{FP.Polynomial{T}}, target))
# end
#
# function GeodesicOnTheSphere{T}(
#     start::MP.AbstractPolynomial,
#     target::MP.AbstractPolynomial) where {T<:Number}
#     s, t = convert(Vector{FP.Polynomial{T}}, [start, target])
#     GeodesicOnTheSphere{T}([s], [t])
# end
#
# function GeodesicOnTheSphere{T}(
#     start::Vector{<:MP.AbstractPolynomial},
#     target::Vector{<:MP.AbstractPolynomial}) where {T<:Number}
#     GeodesicOnTheSphere{T}(
#         convert(Vector{FP.Polynomial{T}}, start),
#         convert(Vector{FP.Polynomial{T}}, target))
# end
#
# function GeodesicOnTheSphere(
#     start::Vector{FP.Polynomial{T}},
#     target::Vector{FP.Polynomial{S}}) where {T<:Number, S<:Number}
#     GeodesicOnTheSphere{promote_type(S,T,Float64)}(start,target)
# end
# function GeodesicOnTheSphere(
#     start::Vector{<:MP.AbstractPolynomial{T}},
#     target::Vector{<:MP.AbstractPolynomial{S}}) where {T<:Number, S<:Number}
#     P = promote_type(S, T, Float64)
#     GeodesicOnTheSphere{P}(
#         convert(Vector{FP.Polynomial{P}}, start),
#         convert(Vector{FP.Polynomial{P}}, target))
# end
# function GeodesicOnTheSphere(
#     start::MP.AbstractPolynomial{T},
#     target::MP.AbstractPolynomial{S}) where {T,S}
#     U = promote_type(S, T, Float64)
#     s, t = convert(Vector{FP.Polynomial{U}}, [start, target])
#     GeodesicOnTheSphere{U}([s], [t])
# end
