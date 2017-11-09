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

const SLH{T} = StraightLineHomotopy{T}

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
function Base.show(io::IO, H::SLH)
    start = join(string.(H.start), ", ")
    target = join(string.(H.target), ", ")
    println(io, typeof(H), "(", "(1-t)⋅[", target, "] + t⋅[", start , "]", ")")
end


function Base.deepcopy(H::SLH)
    StraightLineHomotopy(deepcopy(H.start), deepcopy(H.target))
end

#
# EQUALITY
#
function ==(H1::SLH, H2::SLH)
    H1.start == H2.start && H1.target == H2.target
end
function Base.isequal(H1::SLH, H2::SLH)
    Base.isequal(H1.start, H2.start) && Base.isequal(H1.target, H2.target)
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
function evaluate(H::SLH{T}, x::Vector{S}, t::Number) where {T, S}
    evaluate!(zeros(H.target, promote_type(T, S)), H, x, t)
end
(H::SLH)(x,t) = evaluate(H,x,t)


function evaluate!(u::AbstractVector{T}, H::SLH, x::Vector, t::Number, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}
    evaluate_start_target!(cfg, H, x, precomputed)
    u .= (one(t) - t) .* value_target(cfg) .+ t .* value_start(cfg)
end

function evaluate(H::SLH{T}, x::Vector{S}, t::Number, cfg::PolynomialHomotopyConfig, precomputed=false) where {T, S}
    evaluate!(zeros(H.target, promote_type(T, S)), H, x, t, cfg, precomputed)
end

function jacobian!(u::AbstractMatrix, H::SLH{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}
    jacobian_start_target!(cfg, H, x, precomputed)

    u .= (one(t) - t) .* jacobian_target(cfg) .+ t .* jacobian_start(cfg)
end

function jacobian!(r::JacobianDiffResult, H::SLH{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}
    evaluate_and_jacobian_start_target!(cfg, H, x)

    r.value .= (one(t) - t) .* value_target(cfg) .+ t .* value_start(cfg)
    r.jacobian = (one(t) - t) .* jacobian_target(cfg) .+ t .* jacobian_start(cfg)
    r
end

function jacobian(H::SLH{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}
    u = similar(jacobian_target(cfg))
    jacobian!(u, H, x, t, cfg, precomputed)
    u
end

function dt!(u, H::SLH{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}
    evaluate_start_target!(cfg, H, x, precomputed)
    u .= value_start(cfg) .- value_target(cfg)
end
function dt(H::SLH{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}
    u = similar(value_start(cfg))
    dt!(u, H, x, t, cfg, precomputed)
    u
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
