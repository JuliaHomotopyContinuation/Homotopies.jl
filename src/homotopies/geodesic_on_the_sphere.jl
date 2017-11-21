export GeodesicOnTheSphere

"""
    GeodesicOnTheSphere(start, target)

Homotopy is the geodesic from `g=start/|start|` (t=1) to `f=target/|target|` (t=0):

```H(x,t) = (cos(tα) - sin (tα)cos(α)/sin(α)) f + sin(tα) / sin(α) * g```

where ``α = cos |<f,g>|``. The constructor automatically homgenizes `start` and `target`.

`start` and `target` have to match and to be one of the following
* `Vector{<:MP.AbstractPolynomial}` where `MP` is [`MultivariatePolynomials`](https://github.com/blegat/MultivariatePolynomials.jl)
* `MP.AbstractPolynomial`
* `Vector{<:FP.Polynomial}` where `FP` is [`FixedPolynomials`](https://github.com/saschatimme/FixedPolynomials.jl)


    GeodesicOnTheSphere{T}(start, target)

You can also force a specific coefficient type `T`.
"""
mutable struct GeodesicOnTheSphere{T<:Number} <: AbstractPolynomialHomotopy{T}
    start::Vector{FP.Polynomial{T}}
    target::Vector{FP.Polynomial{T}}
    α::Float64

    # only used for conversions
    function GeodesicOnTheSphere{T}(start::Vector{FP.Polynomial{T}}, target::Vector{FP.Polynomial{T}}, α::Float64) where {T<:Number}
        new(start, target, α)
    end

    function GeodesicOnTheSphere{T}(start::Vector{FP.Polynomial{T}}, target::Vector{FP.Polynomial{T}}) where {T<:Number}
        start, target = homogenize(start, target)

        @assert length(start) == length(target) "Expected the same number of polynomials, but got $(length(start)) and $(length(target))"
        s_nvars = maximum(FP.nvariables.(start))
        @assert all(s_nvars .== FP.nvariables.(start)) "Not all polynomials of the start system have $(s_nvars) variables."
        t_nvars = maximum(FP.nvariables.(target))
        @assert all(t_nvars .== FP.nvariables.(target)) "Not all polynomials of the target system have $(t_nvars) variables instead of $(s_nvars)."
        @assert s_nvars == t_nvars "Expected start and target system to have the same number of variables, but got $(s_nvars) and $(t_nvars)."

        FP.scale_coefficients!.(start, inv(FP.weylnorm(start)))
        FP.scale_coefficients!.(target, inv(FP.weylnorm(target)))

        α = acos(convert(Float64, real(FP.weyldot(start,target))))
        if α > π / 2
            α = π / 2 - α
            FP.scale_coefficients!.(start, -one(T))
        end

        new(start, target, α)
    end

    function GeodesicOnTheSphere{T}(start, target) where {T<:Number}
        s, t = construct(T, start, target)
        GeodesicOnTheSphere{T}(s, t)
    end
end

function GeodesicOnTheSphere(start, target)
    T, s, t = construct(start, target)
    GeodesicOnTheSphere{T}(s, t)
end


function Base.deepcopy(H::GeodesicOnTheSphere{T}) where T
    GeodesicOnTheSphere{T}(deepcopy(H.start), deepcopy(H.target), H.α)
end
#
# PROMOTION AND CONVERSION
#
function Base.promote_rule(::Type{GeodesicOnTheSphere{T}}, ::Type{GeodesicOnTheSphere{S}}) where {S<:Number,T<:Number}
    GeodesicOnTheSphere{promote_type(T,S)}
end
function Base.promote_rule(::Type{GeodesicOnTheSphere{T}}, ::Type{S}) where {S<:Number,T<:Number}
    GeodesicOnTheSphere{promote_type(T,S)}
end
function Base.promote_rule(::Type{GeodesicOnTheSphere},::Type{S}) where {S<:Number}
    GeodesicOnTheSphere{S}
end

function Base.convert(::Type{GeodesicOnTheSphere{T}}, H::GeodesicOnTheSphere) where {T}
    GeodesicOnTheSphere{T}(convert.(FP.Polynomial{T}, H.start), convert.(FP.Polynomial{T}, H.target), H.α)
end


# EVALUATION + DIFFERENTATION
#
(H::GeodesicOnTheSphere)(x,t) = evaluate(H,x,t)

function λ(α, t)
    λ_1 = sin(t * α) / sin(α)
    λ_2 = cos(t * α) - λ_1 * cos(α)

    λ_1, λ_2
end

function dλ(α, t)
    λ_1_dot = α * cos(t * α) / sin(α)
    λ_2_dot = - α * sin(t * α) - λ_1_dot * cos(α)

    λ_1_dot, λ_2_dot
end

function evaluate!(u::AbstractVector, H::GeodesicOnTheSphere{T}, x::Vector, t::Number) where T
    λ_1, λ_2 = λ(H.α, t)
    for i = 1:length(H.target)
        u[i] = λ_1 * FP.evaluate(H.start[i], x) + λ_2 * FP.evaluate(H.target[i], x)
    end
    u
end

function evaluate!(u::AbstractVector{T}, H::GeodesicOnTheSphere, x::Vector, t::Number, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}
    evaluate_start_target!(cfg, H, x, precomputed)
    λ_1, λ_2 = λ(H.α, t)
    u .= λ_1 .* value_start(cfg) .+ λ_2 .* value_target(cfg)
end


function jacobian!(u::AbstractMatrix, H::GeodesicOnTheSphere{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}
    jacobian_start_target!(cfg, H, x, precomputed)
    λ_1, λ_2 = λ(H.α, t)

    u .= λ_1 .* jacobian_start(cfg) .+ λ_2 .* jacobian_target(cfg)
end
function jacobian!(r::JacobianDiffResult, H::GeodesicOnTheSphere{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig) where {T<:Number}
    evaluate_and_jacobian_start_target!(cfg, H, x)
    λ_1, λ_2 = λ(H.α, t)

    r.value .= λ_1 .* value_start(cfg) .+ λ_2 .* value_target(cfg)
    r.jacobian .= λ_1 .* jacobian_start(cfg) .+ λ_2 .* jacobian_target(cfg)
    r
end

function dt!(u, H::GeodesicOnTheSphere{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}
    evaluate_start_target!(cfg, H, x, precomputed)
    λ_1_dot, λ_2_dot = dλ(H.α, t)
    u .= λ_1_dot .* value_start(cfg) .+ λ_2_dot .* value_target(cfg)
end
function dt!(r::DtDiffResult, H::GeodesicOnTheSphere{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}
    evaluate_start_target!(cfg, H, x, precomputed)
    λ_1, λ_2 = λ(H.α, t)
    λ_1_dot, λ_2_dot = dλ(H.α, t)

    r.value = λ_1 .* value_start(cfg) .+ λ_2 .* value_target(cfg)
    r.dt = λ_1_dot .* value_start(cfg) .+ λ_2_dot .* value_target(cfg)
    r
end

function weylnorm(H::GeodesicOnTheSphere{T}) where {T<:Number}
    function (t)
         real(one(T))
     end
end

function gammatrick!(H::GeodesicOnTheSphere{T}, γ::Union{AbstractFloat, Complex}) where {T<:Complex}
    FP.scale_coefficients!.(H.start, convert(T, γ))

    FP.scale_coefficients!.(H.start, inv(FP.weylnorm(H.start)))
    α = acos(convert(Float64, real(FP.weyldot(H.start, H.target))))
    if α > π / 2
        α = π / 2 - α
        FP.scale_coefficients!.(H.start, -one(T))
    end
    H.α = α

    H
end
