export GeodesicOnTheSphere

"""
    GeodesicOnTheSphere(start, target)

Homotopy is the geodesic from g=start/|start| (t=1) to f=target/|target| (t=0):
H(x,t) = (cos(tα) - sin (tα)cos(α)) f + sin(tα) g,
where α = cos <f,g>.

`start` and `target` have to match and to be one of the following
* `Vector{<:MP.AbstractPolynomial}` where `MP` is [`MultivariatePolynomials`](https://github.com/blegat/MultivariatePolynomials.jl)
* `MP.AbstractPolynomial`
* `Vector{<:FP.Polynomial}` where `FP` is [`FixedPolynomials`](https://github.com/saschatimme/FixedPolynomials.jl)


    GeodesicOnTheSphere{T}(start, target)

You can also force a specific coefficient type `T`.
"""
struct GeodesicOnTheSphere{T<:Number} <: AbstractPolynomialHomotopy{T}
    start::Vector{FP.Polynomial{T}}
    target::Vector{FP.Polynomial{T}}

    function GeodesicOnTheSphere{T}(start::Vector{FP.Polynomial{T}}, target::Vector{FP.Polynomial{T}}) where {T<:Number}
        @assert length(start) == length(target) "Expected the same number of polynomials, but got $(length(start)) and $(length(target))"

        s_nvars = maximum(FP.nvariables.(start))
        @assert all(s_nvars .== FP.nvariables.(start)) "Not all polynomials of the start system have $(s_nvars) variables."

        t_nvars = maximum(FP.nvariables.(target))
        @assert all(t_nvars .== FP.nvariables.(target)) "Not all polynomials of the target system have $(t_nvars) variables."

        @assert s_nvars == t_nvars "Expected start and target system to have the same number of variables, but got $(s_nvars) and $(t_nvars)."
        new(start, target)

    end

    function GeodesicOnTheSphere{T}(
        start::Vector{FP.Polynomial{U}},
        target::Vector{FP.Polynomial{V}}
        ) where {T<:Number, U<:Number, V<:Number}
        GeodesicOnTheSphere{T}(
            convert(Vector{FP.Polynomial{T}}, start),
            convert(Vector{FP.Polynomial{T}}, target))
    end

    function GeodesicOnTheSphere{T}(
        start::MP.AbstractPolynomial,
        target::MP.AbstractPolynomial) where {T<:Number}
        s, t = convert(Vector{FP.Polynomial{T}}, [start, target])
        GeodesicOnTheSphere{T}([s], [t])
    end

    function GeodesicOnTheSphere{T}(
        start::Vector{<:MP.AbstractPolynomial},
        target::Vector{<:MP.AbstractPolynomial}) where {T<:Number}
        GeodesicOnTheSphere{T}(
            convert(Vector{FP.Polynomial{T}}, start),
            convert(Vector{FP.Polynomial{T}}, target))
    end
end

#
# CONSTRUCTORS
#
function GeodesicOnTheSphere(
    start::Vector{FP.Polynomial{T}},
    target::Vector{FP.Polynomial{S}}) where {T<:Number, S<:Number}
    GeodesicOnTheSphere{promote_type(S,T)}(start,target)
end
function GeodesicOnTheSphere(
    start::Vector{<:MP.AbstractPolynomial{T}},
    target::Vector{<:MP.AbstractPolynomial{S}}) where {T<:Number, S<:Number}
    P = promote_type(S, T)
    GeodesicOnTheSphere{T}(
        convert(Vector{FP.Polynomial{T}}, start),
        convert(Vector{FP.Polynomial{T}}, target))
end
function GeodesicOnTheSphere(
    start::MP.AbstractPolynomial{T},
    target::MP.AbstractPolynomial{S}) where {T,S}
    U = promote_type(S, T)
    s, t = convert(Vector{FP.Polynomial{T}}, [start, target])
    GeodesicOnTheSphere{T}([s], [t])
end

#
# SHOW
#
function Base.show(io::IO, H::GeodesicOnTheSphere)
    start = join(string.(H.start), ", ")
    target = join(string.(H.target), ", ")
    println(io, typeof(H),
    " The homotopy is given by the spherical geodesic from start/|start| (t=1) to target/|target| (t=0). ",
    "Here: target = [", target, "], start = [", start,"].")
end


#
# EQUALITY
#
function ==(H1::GeodesicOnTheSphere, H2::GeodesicOnTheSphere)
    H1.start == H2.start && H1.target == H2.target
end
function Base.isequal(H1::GeodesicOnTheSphere, H2::GeodesicOnTheSphere)
    Base.isequal(H1.start, H2.start) && Base.isequal(H1.target, H2.target)
end

#
# PROMOTION AND CONVERSION
#
function Base.promote_rule(
    ::Type{GeodesicOnTheSphere{T}},
    ::Type{GeodesicOnTheSphere{S}}) where {S<:Number,T<:Number}
    GeodesicOnTheSphere{promote_type(T,S)}
end

function Base.convert(
    ::Type{GeodesicOnTheSphere{T}},
    H::GeodesicOnTheSphere) where {T}
    GeodesicOnTheSphere{T}(H.start, H.target)
end

#
# EVALUATION + DIFFERENTATION
#
function evaluate!(u::AbstractVector{T}, H::GeodesicOnTheSphere{T}, x::Vector{T}, t::Number) where {T<:Number}
    inverse_norm_start = one(T)/FP.weylnorm(H.start)
    inverse_norm_target = one(T)/FP.weylnorm(H.target)
    α=acos(FP.weyldot(H.start,H.target)  * inverse_norm_start * inverse_norm_target)
    λ_1 = sin(t*α)/sin(α)
    λ_2 = cos(t*α) - λ_1*cos(α)

    map!(u, H.target, H.start) do f, g
        inverse_norm_start * λ_1 * FP.evaluate(g, x) + inverse_norm_target * λ_2 * FP.evaluate(f, x)
    end
end

function evaluate(H::AbstractPolynomialHomotopy{T}, x::Vector{T}, t::Number) where {T<:Number}
    evaluate!(zeros(H.target, T), H, x,  t)
end
(H::GeodesicOnTheSphere)(x,t) = evaluate(H,x,t)

function differentiate(F::Vector{FP.Polynomial{T}}) where {T<:Number}
    [FP.differentiate(f, i) for f in F, i=1:FP.nvariables.(F[1])]
end



function jacobian!(H::GeodesicOnTheSphere{T}) where {T<:Number}
    J_start = differentiate(H.start)
    J_target = differentiate(H.target)
    inverse_norm_start = one(T)/FP.weylnorm(H.start)
    inverse_norm_target = one(T)/FP.weylnorm(H.target)
    α=acos(FP.weyldot(H.start,H.target)  * inverse_norm_start * inverse_norm_target)

    function (u, x, t)
        λ_1 = sin(t*α)/sin(α)
        λ_2 = cos(t*α) - λ_1*cos(α)
        map!(u, J_target, J_start) do f, g
                    inverse_norm_start * λ_1 * FP.evaluate(g, x) + inverse_norm_target * λ_2 * FP.evaluate(f, x)
        end
    end
end
function jacobian(H::GeodesicOnTheSphere{T}) where {T<:Number}
    J_start = differentiate(H.start)
    J_target = differentiate(H.target)
    inverse_norm_start = one(T)/FP.weylnorm(H.start)
    inverse_norm_target = one(T)/FP.weylnorm(H.target)
    α=acos(FP.weyldot(H.start,H.target)  * inverse_norm_start * inverse_norm_target)

    function (x, t)
        λ_1 = sin(t*α)/sin(α)
        λ_2 = cos(t*α) - λ_1*cos(α)
        map(J_target, J_start) do f, g
            inverse_norm_start * λ_1 * FP.evaluate(g, x) + inverse_norm_target * λ_2 * FP.evaluate(f, x)
        end
    end
end

function dt!(H::GeodesicOnTheSphere{T}) where {T<:Number}
    inverse_norm_start = one(T)/FP.weylnorm(H.start)
    inverse_norm_target = one(T)/FP.weylnorm(H.target)
    α=acos(FP.weyldot(H.start,H.target)  * inverse_norm_start * inverse_norm_target)

    function (u, x, t)
        λ_1_dot = t * cos(t*α)/sin(α)
        λ_2_dot = - t * sin(t*α) - λ_1_dot*cos(α)
        map!(u, H.target, H.start) do f, g
            inverse_norm_start * λ_1_dot * FP.evaluate(g, x) + inverse_norm_target * λ_2_dot * FP.evaluate(f, x)
        end
     end
end
function dt(H::GeodesicOnTheSphere{T}) where {T<:Number}
    inverse_norm_start = one(T)/FP.weylnorm(H.start)
    inverse_norm_target = one(T)/FP.weylnorm(H.target)
    α=acos(FP.weyldot(H.start,H.target)  * inverse_norm_start * inverse_norm_target)

    function (x, t)
        λ_1_dot = t * cos(t*α)/sin(α)
        λ_2_dot = - t * sin(t*α) - λ_1_dot*cos(α)
        map(H.target, H.start) do f, g
            inverse_norm_start * λ_1_dot * FP.evaluate(g, x) + inverse_norm_target * λ_2_dot * FP.evaluate(f, x)
        end
     end
end

function homogenize(H::GeodesicOnTheSphere)
    typeof(H)(FP.homogenize.(H.start), FP.homogenize.(H.target))
end
function dehomogenize(H::GeodesicOnTheSphere)
    typeof(H)(FP.dehomogenize.(H.start), FP.dehomogenize.(H.target))
end

function ishomogenized(H::GeodesicOnTheSphere)
    all(FP.ishomogenized.(H.start)) && all(FP.ishomogenized.(H.target))
end
function ishomogenous(H::GeodesicOnTheSphere)
    all(FP.ishomogenous.(H.start)) && all(FP.ishomogenous.(H.target))
end

nvariables(H::GeodesicOnTheSphere) = FP.nvariables(H.start[1])
Base.length(H::GeodesicOnTheSphere) = length(H.start)


#
# """
#     weylnorm(H, t)
#
# Computes the weyl norm of the homotopy `H` to the given time `t`.
#
# ## Explanation
# For ``H = (1-t)F+tG`` we have
# ```math
# \begin{align*}
# <H,H> &= <(1-t)F+tG,(1-t)F+tG> \\
#       &= <(1-t)F,(1-t)F+tG> + <tG,(1-t)F+tG> \\
#       &= <(1-t)F,(1-t)F> + <(1-t)F,tG> + <tG,(1-t)F> + <tG,tG> \\
#       &= <(1-t)F,(1-t)F> + 2real(<(1-t)F,tG>) + <tG,tG> \\
#       &= |1-t|^2<F,F> + 2(t-|t|^2)real(<F,G>) + |t|^2<G,G> \\
# \end{align*}
# ```
# """
# function weylnorm(H::GeodesicOnTheSphere{T}, t::Number) where {T<:Complex}
#     F = H.target
#     G = H.start
#
#     a = abs2(one(T) - t)
#     sqrt(a * FP.weyldot(F,F) + 2 * (t - a) * real(FP.weyldot(F,G)) + abs2(t) * FP.weyldot(G,G))
# end
