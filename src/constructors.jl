function construct(start::Vector{FP.Polynomial{T}}, target::Vector{FP.Polynomial{T}}
    ) where {T<:Number}
    T, start, target
end

function construct(::Type{T}, start::Vector{FP.Polynomial{T}}, target::Vector{FP.Polynomial{T}}
    ) where {T<:Number}
    start, target
end

function construct(::Type{T}, start::Vector{<:FP.Polynomial}, target::Vector{<:FP.Polynomial}
    ) where {T<:Number}
    convert(Vector{FP.Polynomial{T}}, start), convert(Vector{FP.Polynomial{T}}, target)
end


function construct(start::Vector{FP.Polynomial{U}}, target::Vector{FP.Polynomial{V}}
    ) where {U<:Number, V<:Number}
    P = promote_type(U, V, Float64)
    P, convert(Vector{FP.Polynomial{P}}, start), convert(Vector{FP.Polynomial{P}}, target)
end

function construct(::Type{T}, start::Vector{<:MP.AbstractPolynomial}, target::Vector{<:MP.AbstractPolynomial}
    ) where {T<:Number}
    convert(Vector{FP.Polynomial{T}}, start), convert(Vector{FP.Polynomial{T}}, target)
end

function construct(start::Vector{<:MP.AbstractPolynomial{U}}, target::Vector{<:MP.AbstractPolynomial{V}}
    ) where {U<:Number, V<:Number}
    P = promote_type(U, V, Float64)
    P, convert(Vector{FP.Polynomial{P}}, start), convert(Vector{FP.Polynomial{P}}, target)
end

function construct(::Type{T}, start::MP.AbstractPolynomial, target::MP.AbstractPolynomial
    ) where {T<:Number}
    s, t = convert(Vector{FP.Polynomial{T}}, [start, target])
    [s], [t]
end

function construct(start::MP.AbstractPolynomial{U}, target::MP.AbstractPolynomial{V}
    ) where {U<:Number, V<:Number}
    P = promote_type(U, V, Float64)
    s, t = convert(Vector{FP.Polynomial{P}}, [start, target])
    P, [s], [t]
end
