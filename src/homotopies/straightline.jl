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
mutable struct StraightLineHomotopy{T<:Number} <: AbstractPolynomialHomotopy{T}
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

    function StraightLineHomotopy{T}(start, target) where {T<:Number}
        s, t = construct(T, start, target)
        StraightLineHomotopy{T}(s, t)
    end
end

function StraightLineHomotopy(start, target)
    T, s, t = construct(start, target)
    StraightLineHomotopy{T}(s, t)
end

const SLH{T} = StraightLineHomotopy{T}

#
# SHOW
#
function Base.deepcopy(H::SLH)
    StraightLineHomotopy(deepcopy(H.start), deepcopy(H.target))
end

#
# PROMOTION AND CONVERSION
#
Base.promote_rule(::Type{SLH{T}}, ::Type{SLH{S}}) where {S<:Number,T<:Number} = SLH{promote_type(T,S)}
Base.promote_rule(::Type{SLH}, ::Type{S}) where {S<:Number} = SLH{S}
Base.promote_rule(::Type{SLH{T}}, ::Type{S}) where {S<:Number,T<:Number} = SLH{promote_type(T,S)}

Base.convert(::Type{SLH{T}}, H::SLH) where {T} = StraightLineHomotopy{T}(H.start, H.target)


#
# EVALUATION + DIFFERENTATION
#
function evaluate!(u::AbstractVector, H::SLH{T}, x::Vector, t::Number) where T
    for i = 1:length(H.target)
        f = H.target[i]
        g = H.start[i]
        u[i] = (one(T) - t) * FP.evaluate(f, x) + t * FP.evaluate(g, x)
    end
    u
end
(H::SLH)(x,t) = evaluate(H,x,t)


function evaluate!(u::AbstractVector{T}, H::SLH, x::Vector, t::Number, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}
    evaluate_start_target!(cfg, H, x, precomputed)
    u .= (one(t) - t) .* value_target(cfg) .+ t .* value_start(cfg)
end

function jacobian!(u::AbstractMatrix, H::SLH{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}
    jacobian_start_target!(cfg, H, x, precomputed)

    u .= (one(t) - t) .* jacobian_target(cfg) .+ t .* jacobian_start(cfg)
end

function jacobian!(r::JacobianDiffResult, H::SLH{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}
    evaluate_and_jacobian_start_target!(cfg, H, x)

    r.value .= (one(t) - t) .* value_target(cfg) .+ t .* value_start(cfg)
    r.jacobian .= (one(t) - t) .* jacobian_target(cfg) .+ t .* jacobian_start(cfg)
    r
end

function dt!(u, H::SLH{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}
    evaluate_start_target!(cfg, H, x, precomputed)
    u .= value_start(cfg) .- value_target(cfg)
end

function dt!(r::DtDiffResult, H::SLH{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}
    evaluate_start_target!(cfg, H, x, precomputed)
    r.value .= (one(T) - t) .* value_target(cfg) .+ t .* value_start(cfg)
    r.dt .= value_start(cfg) .- value_target(cfg)
    r
end

function weylnorm(H::SLH{T})  where {T<:Number}
    f = FP.homogenize.(H.start)
    g = FP.homogenize.(H.target)
    λ_1 = FP.weyldot(f,f)
    λ_2 = FP.weyldot(f,g)
    λ_3 = FP.weyldot(g,g)

    function (t)
        sqrt(abs2(one(T) - t) * λ_1 + 2 * real((one(T) - t) * conj(t) * λ_2) + abs2(t) * λ_3)
    end
end
