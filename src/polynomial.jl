export gammatrick!

"""
    gammatrick!(H::AbstractPolynomialHomotopy{Complex} [, seed::Int]])

Scale the coefficients of the start system of `H` with a random complex number picked
uniformly from the (complex) unit circle. Use this to make the paths ``z(t)`` generic.

    gammatrick!(H::AbstractPolynomialHomotopy{Complex}, γ::Complex)

You can also pass a scaling factor directly.
"""
function gammatrick!(H::AbstractPolynomialHomotopy{T}, γ::Number) where {T<:Complex}
    FP.scale_coefficients!.(H.start, convert(T, γ))
    H
end
gammatrick!(H::AbstractPolynomialHomotopy{T}) where T = gammatrick!(H, randomgamma(T))
gammatrick!(H::AbstractPolynomialHomotopy{T}, seed::Int) where T = gammatrick!(H, randomgamma(T, seed))

randomgamma(::Type{T}) where {T<:Complex} = convert(T, exp(im * (rand() * 2π - π)))
function randomgamma(::Type{T}, seed::Int) where {T<:Complex}
    srand(seed)
    θ = rand() * 2π - π
    convert(T, exp(im * θ))
end

function homogenize(H::AbstractPolynomialHomotopy)
    typeof(H)(FP.homogenize.(H.start), FP.homogenize.(H.target))
end
function dehomogenize(H::AbstractPolynomialHomotopy)
    typeof(H)(FP.dehomogenize.(H.start), FP.dehomogenize.(H.target))
end

function ishomogenized(H::AbstractPolynomialHomotopy)
    all(FP.ishomogenized.(H.start)) && all(FP.ishomogenized.(H.target))
end
function ishomogenous(H::AbstractPolynomialHomotopy)
    all(FP.ishomogenous.(H.start)) && all(FP.ishomogenous.(H.target))
end

nvariables(H::AbstractPolynomialHomotopy) = FP.nvariables(H.start[1])
Base.length(H::AbstractPolynomialHomotopy) = length(H.start)
