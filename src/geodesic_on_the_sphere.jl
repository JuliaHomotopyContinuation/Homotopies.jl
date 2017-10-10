export GeodesicOnTheSphere

"""
    GeodesicOnTheSphere(start, target)

Homotopy is the geodesic from `g=start/|start|` (t=1) to `f=target/|target|`` (t=0):
``H(x,t) = (cos(tα) - sin (tα)cos(α)) f + sin(tα) g``,
where ``α = cos <f,g>``.

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
    α::Float64

    function GeodesicOnTheSphere{T}(start::Vector{FP.Polynomial{T}}, target::Vector{FP.Polynomial{T}}) where {T<:Number}
        @assert length(start) == length(target) "Expected the same number of polynomials, but got $(length(start)) and $(length(target))"

        s_nvars = maximum(FP.nvariables.(start))
        @assert all(s_nvars .== FP.nvariables.(start)) "Not all polynomials of the start system have $(s_nvars) variables."

        t_nvars = maximum(FP.nvariables.(target))
        @assert all(t_nvars .== FP.nvariables.(target)) "Not all polynomials of the target system have $(t_nvars) variables."

        @assert s_nvars == t_nvars "Expected start and target system to have the same number of variables, but got $(s_nvars) and $(t_nvars)."


        s_norm = FP.weylnorm(FP.homogenize.(start))
        t_norm = FP.weylnorm(FP.homogenize.(target))

        map!(start,start) do f
            FP.Polynomial(f.exponents, f.coefficients./s_norm, f.homogenized)
        end
        map!(target,target) do f
            FP.Polynomial(f.exponents, f.coefficients./t_norm, f.homogenized)
        end

        α = acos(real(FP.weyldot(start,target)))

        if α > π / 2
            α = π / 2 - α
            map!(start, start) do f
                FP.Polynomial(f.exponents, -1.0 .* f.coefficients, f.homogenized)
            end
        end

        new(start, target, α)
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
    GeodesicOnTheSphere{promote_type(S,T,Float64)}(start,target)
end
function GeodesicOnTheSphere(
    start::Vector{<:MP.AbstractPolynomial{T}},
    target::Vector{<:MP.AbstractPolynomial{S}}) where {T<:Number, S<:Number}
    P = promote_type(S, T, Float64)
    GeodesicOnTheSphere{P}(
        convert(Vector{FP.Polynomial{P}}, start),
        convert(Vector{FP.Polynomial{P}}, target))
end
function GeodesicOnTheSphere(
    start::MP.AbstractPolynomial{T},
    target::MP.AbstractPolynomial{S}) where {T,S}
    U = promote_type(S, T, Float64)
    s, t = convert(Vector{FP.Polynomial{U}}, [start, target])
    GeodesicOnTheSphere{U}([s], [t])
end

#
# SHOW
#
function Base.show(io::IO, H::GeodesicOnTheSphere)
    start = join(string.(H.start), ", ")
    target = join(string.(H.target), ", ")
    println(io, typeof(H), "(", "(1-t)⋅[", target, "] + t⋅[", start , "]", ") with angle $(H.α)")
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
    λ_1 = sin(t*H.α)/sin(H.α)
    λ_2 = cos(t*H.α) - λ_1*cos(H.α)
    map!(u, H.target, H.start) do f, g
        λ_1 * FP.evaluate(g, x) + λ_2 * FP.evaluate(f, x)
    end
end

function evaluate(H::AbstractPolynomialHomotopy{T}, x::Vector{T}, t::Number) where {T<:Number}
    evaluate!(zeros(H.target, T), H, x,  t)
end
(H::GeodesicOnTheSphere)(x,t) = evaluate(H,x,t)


function weylnorm(H::GeodesicOnTheSphere{Complex{T}})  where {T<:Real}
    function (t)
         one(T)
     end
end


function differentiate(F::Vector{FP.Polynomial{T}}) where {T<:Number}
    [FP.differentiate(f, i) for f in F, i=1:FP.nvariables.(F[1])]
end



function jacobian!(H::GeodesicOnTheSphere{T}) where {T<:Number}
    J_start = differentiate(H.start)
    J_target = differentiate(H.target)

    function (u, x, t)
        λ_1 = sin(t * H.α) / sin(H.α)
        λ_2 = cos(t * H.α) - λ_1 * cos(H.α)
        map!(u, J_target, J_start) do f, g
            λ_1 * FP.evaluate(g, x) +  λ_2 * FP.evaluate(f, x)
        end
    end
end

function jacobian(H::GeodesicOnTheSphere{T}) where {T<:Number}
    J_start = differentiate(H.start)
    J_target = differentiate(H.target)

    function (x, t)
        λ_1 = sin(t * H.α) / sin(H.α)
        λ_2 = cos(t * H.α) - λ_1 * cos(H.α)
        map(J_target, J_start) do f, g
             λ_1 * FP.evaluate(g, x) + λ_2 * FP.evaluate(f, x)
        end
    end
end

function dt!(H::GeodesicOnTheSphere{T}) where {T<:Number}

    function (u, x, t)
        λ_1_dot = t * cos(t * H.α)/sin(H.α)
        λ_2_dot = - t * sin(t * H.α) - λ_1_dot * cos(H.α)
        map!(u, H.target, H.start) do f, g
            λ_1_dot * FP.evaluate(g, x) + λ_2_dot * FP.evaluate(f, x)
        end
     end
end
function dt(H::GeodesicOnTheSphere{T}) where {T<:Number}

    function (x, t)
        λ_1_dot = t * cos(t * H.α) / sin(H.α)
        λ_2_dot = - t * sin(t * H.α) - λ_1_dot * cos(H.α)
        map(H.target, H.start) do f, g
            λ_1_dot * FP.evaluate(g, x) + λ_2_dot * FP.evaluate(f, x)
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
