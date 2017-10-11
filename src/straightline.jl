export StraightLineHomotopy

"""
    StraightLineHomotopy(start, target)

Construct the homotopy `t * start + (1-t) * target`.

`start` and `target` have to match and to be one of the following
* `Vector{<:MP.AbstractPolynomial}` where `MP` is [`MultivariatePolynomials`](https://github.com/blegat/MultivariatePolynomials.jl)
* `MP.AbstractPolynomial`
* `Vector{<:FP.Polynomial}` where `FP` is [`FixedPolynomials`](https://github.com/saschatimme/FixedPolynomials.jl)


    StraightLineHomotopy{T}(start, target)

You can also force a specific coefficient type `T`.
"""
struct StraightLineHomotopy{T<:Number} <: AbstractPolynomialHomotopy{T}
    start::Vector{FP.Polynomial{T}}
    target::Vector{FP.Polynomial{T}}

    function StraightLineHomotopy{T}(start::Vector{FP.Polynomial{T}}, target::Vector{FP.Polynomial{T}}) where {T<:Number}
        @assert length(start) == length(target) "Expected the same number of polynomials, but got $(length(start)) and $(length(target))"

        s_nvars = maximum(FP.nvariables.(start))
        @assert all(s_nvars .== FP.nvariables.(start)) "Not all polynomials of the start system have $(s_nvars) variables."

        t_nvars = maximum(FP.nvariables.(target))
        @assert all(t_nvars .== FP.nvariables.(target)) "Not all polynomials of the target system have $(t_nvars) variables."

        @assert s_nvars == t_nvars "Expected start and target system to have the same number of variables, but got $(s_nvars) and $(t_nvars)."
        new(start, target)
    end

    function StraightLineHomotopy{T}(
        start::Vector{FP.Polynomial{U}},
        target::Vector{FP.Polynomial{V}}
        ) where {T<:Number, U<:Number, V<:Number}
        StraightLineHomotopy{T}(
            convert(Vector{FP.Polynomial{T}}, start),
            convert(Vector{FP.Polynomial{T}}, target))
    end

    function StraightLineHomotopy{T}(
        start::MP.AbstractPolynomial,
        target::MP.AbstractPolynomial) where {T<:Number}
        s, t = convert(Vector{FP.Polynomial{T}}, [start, target])
        StraightLineHomotopy{T}([s], [t])
    end

    function StraightLineHomotopy{T}(
        start::Vector{<:MP.AbstractPolynomial},
        target::Vector{<:MP.AbstractPolynomial}) where {T<:Number}
        StraightLineHomotopy{T}(
            convert(Vector{FP.Polynomial{T}}, start),
            convert(Vector{FP.Polynomial{T}}, target))
    end
end

#
# CONSTRUCTORS
#
function StraightLineHomotopy(
    start::Vector{FP.Polynomial{T}},
    target::Vector{FP.Polynomial{S}}) where {T<:Number, S<:Number}
    StraightLineHomotopy{promote_type(S,T)}(start,target)
end
function StraightLineHomotopy(
    start::Vector{<:MP.AbstractPolynomial{T}},
    target::Vector{<:MP.AbstractPolynomial{S}}) where {T<:Number, S<:Number}
    P = promote_type(S, T)
    StraightLineHomotopy{P}(
        convert(Vector{FP.Polynomial{T}}, start),
        convert(Vector{FP.Polynomial{T}}, target))
end
function StraightLineHomotopy(
    start::MP.AbstractPolynomial{T},
    target::MP.AbstractPolynomial{S}) where {T,S}
    U = promote_type(S, T)
    s, t = convert(Vector{FP.Polynomial{U}}, [start, target])
    StraightLineHomotopy{U}([s], [t])
end

#
# SHOW
#
function Base.show(io::IO, H::StraightLineHomotopy)
    start = join(string.(H.start), ", ")
    target = join(string.(H.target), ", ")
    println(io, typeof(H), "(", "(1-t)⋅[", target, "] + t⋅[", start , "]", ")")
end

#
# EQUALITY
#
function ==(H1::StraightLineHomotopy, H2::StraightLineHomotopy)
    H1.start == H2.start && H1.target == H2.target
end
function Base.isequal(H1::StraightLineHomotopy, H2::StraightLineHomotopy)
    Base.isequal(H1.start, H2.start) && Base.isequal(H1.target, H2.target)
end

#
# PROMOTION AND CONVERSION
#
function Base.promote_rule(
    ::Type{StraightLineHomotopy{T}},
    ::Type{StraightLineHomotopy{S}}) where {S<:Number,T<:Number}
    StraightLineHomotopy{promote_type(T,S)}
end

function Base.convert(
    ::Type{StraightLineHomotopy{T}},
    H::StraightLineHomotopy) where {T}
    StraightLineHomotopy{T}(H.start, H.target)
end

#
# EVALUATION + DIFFERENTATION
#
function evaluate!(u::AbstractVector{T}, H::StraightLineHomotopy{T}, x::Vector{T}, t::Number) where {T<:Number}
    map!(u, H.target, H.start) do f, g
        (one(T) - t) * FP.evaluate(f, x) + t * FP.evaluate(g, x)
    end
end
function evaluate(H::StraightLineHomotopy{T}, x::Vector{T}, t::Number) where {T<:Number}
    evaluate!(zeros(H.target, T), H, x,  t)
end
(H::StraightLineHomotopy)(x,t) = evaluate(H,x,t)


function weylnorm(H::StraightLineHomotopy{T})  where {T<:Number}
    f = FP.homogenize.(H.start)
    g = FP.homogenize.(H.target)
    λ_1 = FP.weyldot(f,f)
    λ_2 = FP.weyldot(f,g)
    λ_3 = FP.weyldot(g,g)

    function (t)
        sqrt(abs2(one(T) - t) * λ_1 + 2 * real((one(T) - t) * conj(t) * λ_2) + abs2(t) * λ_3)
    end
end

function jacobian!(H::StraightLineHomotopy{T}) where {T<:Number}
    J_start = differentiate(H.start)
    J_target = differentiate(H.target)

    function (u, x, t)
        map!(u, J_target, J_start) do f, g
            (one(T) - t) * FP.evaluate(f, x) + t * FP.evaluate(g, x)
        end
    end
end
function jacobian(H::StraightLineHomotopy{T}) where {T<:Number}
    J_start = differentiate(H.start)
    J_target = differentiate(H.target)

    function (x, t)
        map(J_target, J_start) do f, g
            (one(T) - t) * FP.evaluate(f, x) + t * FP.evaluate(g, x)
        end
    end
end

function dt!(H::StraightLineHomotopy{T}) where {T<:Number}
    function (u, x, ::Number)
        map!(u, H.target, H.start) do f, g
            FP.evaluate(g, x) - FP.evaluate(f, x)
        end
     end
end
function dt(H::StraightLineHomotopy{T}) where {T<:Number}
    function (x, ::Number)
        map(H.target, H.start) do f, g
            FP.evaluate(g, x) - FP.evaluate(f, x)
        end
     end
end

function homogenize(H::StraightLineHomotopy)
    typeof(H)(FP.homogenize.(H.start), FP.homogenize.(H.target))
end
function dehomogenize(H::StraightLineHomotopy)
    typeof(H)(FP.dehomogenize.(H.start), FP.dehomogenize.(H.target))
end

function ishomogenized(H::StraightLineHomotopy)
    all(FP.ishomogenized.(H.start)) && all(FP.ishomogenized.(H.target))
end
function ishomogenous(H::StraightLineHomotopy)
    all(FP.ishomogenous.(H.start)) && all(FP.ishomogenous.(H.target))
end

nvariables(H::StraightLineHomotopy) = FP.nvariables(H.start[1])
Base.length(H::StraightLineHomotopy) = length(H.start)
