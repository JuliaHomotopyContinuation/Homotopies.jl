export GammaTrickHomotopy

"""
    GammaTrickHomotopy(start, target, [γ])

Construct the homotopy `(1-t)⋅target + t⋅γ⋅start`. If `γ` is not supplied it will be drawn
randomly and uniformly from the complex unit circle.

`start` and `target` have to match and to be one of the following
* `Vector{<:MP.AbstractPolynomial}` where `MP` is [`MultivariatePolynomials`](https://github.com/blegat/MultivariatePolynomials.jl)
* `MP.AbstractPolynomial`
* `Vector{<:FP.Polynomial}` where `FP` is [`FixedPolynomials`](https://github.com/saschatimme/FixedPolynomials.jl)


    GammaTrickHomotopy(start, target, seed::Int)

You can also supply a `seed` for the RNG which is used to draw `γ`,
i.e. subsequent invocations with the same `seed` will produce the same `γ`.


    GammaTrickHomotopy{T}(start, target, [γ])
    GammaTrickHomotopy{T}(start, target, seed::Int)

You can also force a specific coefficient type `T`.
"""
struct GammaTrickHomotopy{T<:Complex} <: AbstractPolynomialHomotopy{T}
    start::Vector{FP.Polynomial{T}}
    target::Vector{FP.Polynomial{T}}
    γ::T


    function GammaTrickHomotopy{T}(start::Vector{FP.Polynomial{T}}, target::Vector{FP.Polynomial{T}}, γ::T) where {T<:Complex}
        @assert length(start) == length(target) "Expected the same number of polynomials, but got $(length(start)) and $(length(target))"

        s_nvars = maximum(FP.nvariables.(start))
        @assert all(s_nvars .== FP.nvariables.(start)) "Not all polynomials of the start system have $(s_nvars) variables."

        t_nvars = maximum(FP.nvariables.(target))
        @assert all(t_nvars .== FP.nvariables.(target)) "Not all polynomials of the target system have $(t_nvars) variables."

        @assert s_nvars == t_nvars "Expected start and target system to have the same number of variables, but got $(s_nvars) and $(t_nvars)."
        new(start, target, γ)
    end

    function GammaTrickHomotopy{T}(
        start::Vector{FP.Polynomial{U}},
        target::Vector{FP.Polynomial{V}},
        γ::Complex=randomgamma(T)) where {T<:Complex, U<:Number, V<:Number}
        GammaTrickHomotopy{T}(
            convert(Vector{FP.Polynomial{T}}, start),
            convert(Vector{FP.Polynomial{T}}, target),
            convert(T, γ))
    end

    function GammaTrickHomotopy{T}(
        start::MP.AbstractPolynomial,
        target::MP.AbstractPolynomial,
        γ::Complex=randomgamma(T)) where {T<:Complex}
        s, t = convert(Vector{FP.Polynomial{T}}, [start, target])
        GammaTrickHomotopy{T}([s], [t], convert(T, γ))
    end

    function GammaTrickHomotopy{T}(
        start::Vector{<:MP.AbstractPolynomial{U}},
        target::Vector{<:MP.AbstractPolynomial{V}},
        γ::Complex=randomgamma(T)) where {T<:Complex, U<:Number, V<:Number}
        GammaTrickHomotopy{T}(
            convert(Vector{FP.Polynomial{T}}, start),
            convert(Vector{FP.Polynomial{T}}, target),
            convert(T, γ))
    end

    function GammaTrickHomotopy{T}(start, target, seed::Int) where {T<:Complex}
        GammaTrickHomotopy{T}(start, target, randomgamma(T, seed))
    end
end

#
# CONSTRUCTORS
#
function GammaTrickHomotopy(
    start::Vector{FP.Polynomial{T}},
    target::Vector{FP.Polynomial{S}},
    γ=randomgamma(_promote_type(S,T))) where {T<:Number, S<:Number}
    U = _promote_type(S,T)
    GammaTrickHomotopy{U}(start, target, convert(U, γ))
end

function GammaTrickHomotopy(
    start::Vector{FP.Polynomial{T}},
    target::Vector{FP.Polynomial{S}},
    seed::Int) where {T<:Number, S<:Number}
    U = _promote_type(S,T)
    GammaTrickHomotopy{U}(start, target, randomgamma(U, seed))
end

function GammaTrickHomotopy(
    start::Vector{<:MP.AbstractPolynomial{T}},
    target::Vector{<:MP.AbstractPolynomial{S}},
    γ=randomgamma(_promote_type(S, T))) where {T<:Number, S<:Number}
    P = _promote_type(S, T)
    GammaTrickHomotopy{P}(convert(Vector{FP.Polynomial{P}}, start), convert(Vector{FP.Polynomial{P}}, target), convert(P, γ))
end

function GammaTrickHomotopy(
    start::Vector{<:MP.AbstractPolynomial{T}},
    target::Vector{<:MP.AbstractPolynomial{S}},
    seed::Int) where {T<:Complex, S<:Number}
    P = _promote_type(S, T)
    GammaTrickHomotopy{P}(convert(Vector{FP.Polynomial{P}}, start), convert(Vector{FP.Polynomial{P}}, target), randomgamma(P, seed))
end

function GammaTrickHomotopy(start::MP.AbstractPolynomial{T}, target::MP.AbstractPolynomial{S}) where {T<:Number,S<:Number}
    U = _promote_type(S, T)
    s, t = convert(Vector{FP.Polynomial{U}}, [start, target])
    GammaTrickHomotopy{U}([s], [t], randomgamma(U))
end


_promote_type(::Type{T}, ::Type{S}) where {T<:Number, S<:Number} = complex(float(promote_type(S, T)))

randomgamma(::Type{T}) where {T<:Complex} = convert(T, exp(im * (rand() * 2π - π)))
function randomgamma(::Type{T}, seed::Int) where {T<:Complex}
    srand(seed)
    θ = rand() * 2π - π
    convert(T, exp(im * θ))
end

#
# SHOW
#
function Base.show(io::IO, H::GammaTrickHomotopy)
    start = join(string.(H.start), ", ")
    target = join(string.(H.target), ", ")
    println(io, typeof(H), "(", "(1-t)⋅[", target, "] + t⋅[", start, "]", ")")
end

#
# EQUALITY
#
function ==(H1::GammaTrickHomotopy, H2::GammaTrickHomotopy)
    H1.start == H2.start && H1.target == H2.target && H1.γ == H2.γ
end
function Base.isequal(H1::GammaTrickHomotopy, H2::GammaTrickHomotopy)
    Base.isequal(H1.start, H2.start) &&
    Base.isequal(H1.target, H2.target) &&
    Base.isequal(H1.γ, H2.γ)
end

#
# PROMOTION AND CONVERSION
#
function Base.promote_rule(
    ::Type{GammaTrickHomotopy{T}},
    ::Type{GammaTrickHomotopy{S}}) where {S<:Complex,T<:Complex}
    GammaTrickHomotopy{promote_type(T,S)}
end

function Base.convert(
    ::Type{GammaTrickHomotopy{T}},
    H::GammaTrickHomotopy) where {T}
    GammaTrickHomotopy{T}(H.start, H.target, convert(T, H.γ))
end

#
# EVALUATION + DIFFERENTATION
#
function evaluate!(u::AbstractVector{T}, H::GammaTrickHomotopy{T}, x::Vector{T}, t::Number) where {T<:Complex}
    map!(u, H.target, H.start) do f, g
        (one(T) - t) * FP.evaluate(f, x) + t * H.γ * FP.evaluate(g, x)
    end
end
function evaluate(H::AbstractPolynomialHomotopy{T}, x::Vector{T}, t::Number) where {T<:Complex}
    evaluate!(zeros(H.target, T), H, x,  t)
end
(H::GammaTrickHomotopy)(x,t) = evaluate(H,x,t)

function differentiate(F::Vector{FP.Polynomial{T}}) where {T<:Complex}
    [FP.differentiate(f, i) for f in F, i=1:FP.nvariables.(F[1])]
end
function jacobian!(H::GammaTrickHomotopy{T}) where {T<:Complex}
    J_start = differentiate(H.start)
    J_target = differentiate(H.target)
    γ = H.γ

    function (u, x, t::Number)
        map!(u, J_target, J_start) do f, g
            (1 - t) * FP.evaluate(f, x) + γ * t * FP.evaluate(g, x)
        end
    end
end
function jacobian(H::GammaTrickHomotopy{T}) where {T<:Complex}
    J_start = differentiate(H.start)
    J_target = differentiate(H.target)
    γ = H.γ

    function (x, t::Number)
        map(J_target, J_start) do f, g
            (1 - t) * FP.evaluate(f, x) + t * γ * FP.evaluate(g, x)
        end
    end
end

function dt!(H::GammaTrickHomotopy{T}) where {T<:Complex}
    γ = H.γ
    function (u, x, ::Number)
        map!(u, H.target, H.start) do f, g
            γ * FP.evaluate(g, x) - FP.evaluate(f, x)
        end
     end
end
function dt(H::GammaTrickHomotopy{T}) where {T<:Complex}
    γ = H.γ
    function (x, ::Number)
        map(H.target, H.start) do f, g
            γ * FP.evaluate(g, x) - FP.evaluate(f, x)
        end
     end
end


function homogenize(H::GammaTrickHomotopy)
    typeof(H)(FP.homogenize.(H.start), FP.homogenize.(H.target), H.γ)
end
function dehomogenize(H::GammaTrickHomotopy)
    typeof(H)(FP.dehomogenize.(H.start), FP.dehomogenize.(H.target), H.γ)
end

function ishomogenized(H::GammaTrickHomotopy)
    all(FP.ishomogenized.(H.start)) && all(FP.ishomogenized.(H.target))
end
function ishomogenous(H::GammaTrickHomotopy)
    all(FP.ishomogenous.(H.start)) && all(FP.ishomogenous.(H.target))
end

nvariables(H::GammaTrickHomotopy) = FP.nvariables(H.start[1])
Base.length(H::GammaTrickHomotopy) = length(H.start)

# """
#     weylnorm(H, t)
#
# Computes the weyl norm of the homotopy `H` to the given time `t`.
#
# ## Explanation
# For ``H = (1-t)F+tG`` we have
# ```math
# \begin{align*}
# <H,H> &= <(1-t)F+tγG,(1-t)F+tγG> \\
#       &= <(1-t)F,(1-t)F+tγG> + <tγG,(1-t)F+tγG> \\
#       &= <(1-t)F,(1-t)F> + <(1-t)F,tγG> + <tγG,(1-t)F> + <tγG,tγG> \\
#       &= <(1-t)F,(1-t)F> + 2real(<(1-t)F,tγG>) + <tγG,tγG> \\
#       &= |1-t|^2<F,F> + 2real(γ(t-|t|^2)<F,G>) + |γ|^2|t|^2<G,G>
# \end{align*}
# ```
# """
# function weylnorm(H::GammaTrickHomotopy{T}, t::Number) where {T<:Complex}
#     F = H.target
#     G = H.start
#
#     a = abs2(one(T) - t)
#
#     sqrt(a * FP.weyldot(F,F) + 2 * real(H.γ * (t - a) * FP.weyldot(F,G)) + abs2(H.γ) * abs2(t) * FP.weyldot(G,G))
# end
