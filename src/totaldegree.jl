export TotalDegreeSolutionIterator, totaldegree_startsystem, totaldegree

"""
    TotalDegreeSolutionIterator(degrees, b)

Given the `Vector`s `degrees` and `b` `TotalDegreeSolutionIterator` enumerates all solutions
of the system
```math
\\begin{align*}
    z_1^{d_1} - b_1 &= 0 \\\\
    z_1^{d_2} - b_2 &= 0 \\\\
    &\\vdots \\\\
    z_n^{d_n} - b_n &= 0 \\\\
\\end{align*}
```
where ``d_i`` is `degrees[i]` and ``b_i`` is `b[i]`.
"""
struct TotalDegreeSolutionIterator{T<:Complex, Iter}
    degrees::Vector{Int}
    b::Vector{T}

    iterator::Iter
end
function TotalDegreeSolutionIterator(degrees::Vector{Int}, b::Vector{T}) where {T<:Complex}
    iterator = Base.Iterators.product(map(d -> 0:d-1, degrees)...)
    TotalDegreeSolutionIterator(degrees, b, iterator)
end

Base.start(iter::TotalDegreeSolutionIterator) = start(iter.iterator)
function Base.next(iter::TotalDegreeSolutionIterator{T}, state) where {T<:Complex}
    indices, nextstate = next(iter.iterator, state)

    value = map(1:length(indices), indices) do i, k
        d = iter.degrees[i]
        convert(T, iter.b[i]^(1/d) * exp(2π*im/d)^k)
    end
    value, nextstate
end
Base.done(iter::TotalDegreeSolutionIterator, state) = done(iter.iterator, state)
Base.length(iter::TotalDegreeSolutionIterator) = length(iter.iterator)
Base.eltype(iter::TotalDegreeSolutionIterator{T}) where {T} = Vector{T}

"""
    totaldegree(H::Type{AbstractPolynomialHomotopy}, F, [unitroots=false])

Construct a  total degree homotopy of type `H` with `F` and an iterator of its solutions.
This is the homotopy with start system
```math
\\begin{align*}
    z_1^{d_1} &- b_1\\\\
    z_1^{d_2} &- b_2\\\\
    &\\vdots \\\\
    z_n^{d_n} &- b_n\\\\
\\end{align*}
```
and target system `F`, where ``d_i`` is the degree of the ``i``-th polynomial of `F`.
If `unitroots=true` then ``b_i=1`` otherwise ``b_i`` is a random
complex number (with real and imaginary part in the unit interval).

## Example
```julia
H, startsolutions = totaldegree(StraightLineHomotopy{Complex128}, [x^2+y+1, x^3*y-2])
```
"""
function totaldegree(H::Type{<:AbstractPolynomialHomotopy{T}}, fs::Vector{<:MP.AbstractPolynomial{U}}; kwargs...) where {T<:Complex, U<:Number}
    F = convert(Vector{FP.Polynomial{T}}, fs)
    totaldegree(H, F; kwargs...)
end
function totaldegree(H::Type{<:AbstractPolynomialHomotopy}, fs::Vector{<:MP.AbstractPolynomial{U}}; kwargs...) where {U<:Number}
    F = convert(Vector{FP.Polynomial{promote_type(Complex128, U)}}, fs)
    totaldegree(H, F; kwargs...)
end

function totaldegree(H::Type{<:AbstractPolynomialHomotopy}, F::Vector{FP.Polynomial{U}}; kwargs...) where {U<:Number}
    G, solutions = totaldegree_startsystem(F; kwargs...)
    H(G, F), solutions
end


"""
    totaldegree_startsystem(F::Vector{FP.Polynomial{<:Complex}}, [unit_roots=false])

Return the system
```math
\\begin{align*}
    z_1^{d_1} &- b_1\\\\
    z_1^{d_2} &- b_2\\\\
    &\\vdots \\\\
    z_n^{d_n} &- b_n\\\\
\\end{align*}
```
where ``d_i`` is the degree of the ``i``-th polynomial of `F` and an iterator of its
solutions.
If `unitroots=true` then ``b_i=1`` otherwise ``b_i`` is a random
complex number (with real and imaginary part in the unit interval).
"""
function totaldegree_startsystem(F::Vector{FP.Polynomial{Complex{T}}}; unitroots=false) where {T<:AbstractFloat}
    vars = FP.variables(F[1])
    degrees = FP.degree.(F)
    if unitroots
        b = ones(Complex{T}, length(degrees))
    else
        b = convert.(Complex{T}, rand(Complex128, length(degrees)))
    end

    n = length(degrees)
    N = length(vars)
    make_homogenous = !all(FP.ishomogenous, F) && length(vars) == length(F) + 1
    if n + 1 == N
        if !all(FP.ishomogenous, F)
            return error("In order to create a total degree start system your input system needs " *
                "to have the same number of variables as number of equations or to have N variables
                and N - 1 homogenous polynomials. Currently your system has " *
                "$(length(vars)) variables and $(length(degrees)) equations and is not homogenous.")
        end

        G = map(2:N, degrees) do i, d
            exponents = zeros(Int, N, 2)
            exponents[i, 1] = d
            exponents[1, 2] = d
            coeffs = zeros(Complex{T}, 2)
            coeffs[1] = one(Complex{T})
            coeffs[2] = -b[i - 1]
            FP.Polynomial(exponents, coeffs, vars)
        end

        startsolutions::Vector{Vector{Complex{T}}} =
            map(TotalDegreeSolutionIterator(degrees, b)) do sol
                unshift!(sol, 1)
            end
        G, startsolutions
    elseif n != N
        return error("In order to create a total degree start system your input system needs " *
            "to have the same number of variables as number of equations. Currently your system has " *
            "$(length(vars)) variables and $(length(degrees)) equations.")
    else
        G = map(1:n, degrees) do i, d
            exponents = zeros(Int, n, 2)
            exponents[i, 1] = d
            coeffs = zeros(Complex{T}, 2)
            coeffs[1] = one(Complex{T})
            coeffs[2] = -b[i]
            FP.Polynomial(exponents, coeffs, vars)
        end
        G, collect(TotalDegreeSolutionIterator(degrees, b))::Vector{Vector{Complex{T}}}
    end
end
function totaldegree_startsystem(F::Vector{FP.Polynomial{T}}; kwargs...) where {T<:Number}
    totaldegree_startsystem(convert.(FP.Polynomial{complex(float(T))}, F); kwargs...)
end
