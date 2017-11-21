export gammatrick!, gammatrick


#
# EVALUATION + DIFFERENTATION
#
function evaluate(H::AbstractPolynomialHomotopy{T}, x::Vector{S}, t::Number) where {T, S}
    evaluate!(zeros(H.target, promote_type(T, S)), H, x, t)
end


function evaluate(H::AbstractPolynomialHomotopy{T}, x::Vector{S}, t::Number, cfg::PolynomialHomotopyConfig, precomputed=false) where {T, S}
    evaluate!(zeros(H.target, promote_type(T, S)), H, x, t, cfg, precomputed)
end


function jacobian(H::AbstractPolynomialHomotopy{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}
    u = similar(jacobian_target(cfg))
    jacobian!(u, H, x, t, cfg, precomputed)
    u
end

function dt(H::AbstractPolynomialHomotopy{T}, x::AbstractVector, t, cfg::PolynomialHomotopyConfig, precomputed=false) where {T<:Number}
    u = similar(value_start(cfg))
    dt!(u, H, x, t, cfg, precomputed)
    u
end

#
# EQUALITY
#
function ==(H1::T, H2::T) where {T<:AbstractPolynomialHomotopy}
    H1.start == H2.start && H1.target == H2.target
end
function Base.isequal(H1::T, H2::T) where {T<:AbstractPolynomialHomotopy}
    Base.isequal(H1.start, H2.start) && Base.isequal(H1.target, H2.target)
end

function Base.show(io::IO, H::AbstractPolynomialHomotopy)
    println(io, typeof(H), " with $(length(H.start)) polynomials in $(nvariables(H)) variables.")
end



"""
    gammatrick!(H::AbstractPolynomialHomotopy{Complex} [, seed::Int]])

Scale the coefficients of the start system of `H` with a random complex number picked
uniformly from the (complex) unit circle. Use this to make the paths ``z(t)`` generic.

    gammatrick!(H::AbstractPolynomialHomotopy{Complex}, γ::Complex)

You can also pass a scaling factor directly.
"""
function gammatrick!(H::AbstractPolynomialHomotopy{T}, γ::Union{AbstractFloat, Complex}) where {T<:Complex}
    FP.scale_coefficients!.(H.start, convert(T, γ))
    γ
end
gammatrick!(H::AbstractPolynomialHomotopy{T}) where T = gammatrick!(H, randomgamma(T))
gammatrick!(H::AbstractPolynomialHomotopy{T}, seed::Int) where T = gammatrick!(H, randomgamma(T, seed))

"""
    gammatrick(H::AbstractPolynomialHomotopy{Complex} , γ::Number)

Scale the coefficients of the start system of `H` with `γ`.

    gammatrick(H::AbstractPolynomialHomotopy{Complex})

A a random complex number `γ` is picked uniformly from the (complex) unit circle and then
scale the coefficients of the start system of `H` with `γ`.
This returns the new `H` and `γ`.
"""
function gammatrick(H::AbstractPolynomialHomotopy)
    K = deepcopy(H)
    γ = gammatrick!(K)
    K, γ
end
function gammatrick(H::AbstractPolynomialHomotopy, γ::Number)
    K = deepcopy(H)
    gammatrick!(K, γ)
    K
end

randomgamma(::Type{T}) where {T<:Complex} = convert(T, exp(im * (rand() * 2π - π)))
function randomgamma(::Type{T}, seed::Int) where {T<:Complex}
    srand(seed)
    θ = rand() * 2π - π
    convert(T, exp(im * θ))
end


function homogenize(F::Vector{<:FP.Polynomial}, G::Vector{<:FP.Polynomial})
    F_ishomogenous = all(FP.ishomogenous, F)
    G_ishomogenous = all(FP.ishomogenous, G)

    if F_ishomogenous && G_ishomogenous
        F1 = map(f -> FP.homogenize(f, respect_homogenous=true), F)
        G1 = map(g -> FP.homogenize(g, respect_homogenous=true), G)
        (F1, G1)
    else
        F1 = map(f -> FP.homogenize(f, respect_homogenous=false), F)
        G1 = map(g -> FP.homogenize(g, respect_homogenous=false), G)
        (F1, G1)
    end
end

function homogenize(H::AbstractPolynomialHomotopy)
    G, F = homogenize(H.start, H.target)
    typeof(H)(G, F)
end
function dehomogenize(H::AbstractPolynomialHomotopy)
    typeof(H)(FP.dehomogenize.(H.start), FP.dehomogenize.(H.target))
end

function ishomogenized(H::AbstractPolynomialHomotopy)
    all(FP.ishomogenized, H.start) && all(FP.ishomogenized, H.target)
end
function ishomogenous(H::AbstractPolynomialHomotopy)
    all(FP.ishomogenous, H.start) && all(FP.ishomogenous, H.target)
end

nvariables(H::AbstractPolynomialHomotopy) = FP.nvariables(H.start[1])
Base.length(H::AbstractPolynomialHomotopy) = length(H.start)
