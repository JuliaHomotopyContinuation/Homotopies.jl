export randomhomotopy, randomsystem

"""
    randomhomotopy(::Type{AbstractPolynomialHomotopy{T}}, size::Int; kwargs...)

Create a total degree homotopy where the target system is a [`randomsystem(T, size, size; kwargs...)`](@ref).

## Example
```julia-repl
julia> H, solutions = randomhomotopy(StraightLineHomotopy{Complex128}, 2, mindegree=3, maxdegree=6);
julia> length(H)
3
julia> nvariables(H)
3
```
"""
function randomhomotopy(H::Type{<:AbstractPolynomialHomotopy{T}}, size::Int; kwargs...) where {T<:Complex}
    F = randomsystem(T, size, size; kwargs...)
    totaldegree(H, F)
end
function randomhomotopy(H::Type{<:AbstractPolynomialHomotopy}, size::Int; kwargs...)
    F = randomsystem(size, size; kwargs...)
    totaldegree(H, F)
end

"""
    randomsystem([T=Complex128,] nequations::Int, nvars::Int; mindegree=0, maxdegree=5, rng=Base.Random.GLOBAL_RNG)

Creates a random polynomial system of `nequations` equations with `nvars` variables (named ``x_1``, ...``x_{nvars}``).
Each polynomial has a total degree uniformly drawn from ``[mindegree, maxdegree]``.
The coefficients are drawn independently from the given `rng`.

    randomsystem([T=Complex128,] degrees::Vector{Int}, variables::Vector{Symbol}; rng=N(0,1))

Create a random polynomial system with the given `degrees` and `variables`.
"""
function randomsystem(::Type{T}, nequations::Int, nvars::Int; mindegree=2, maxdegree=5, rng=Base.Random.GLOBAL_RNG) where {T<:Number}
    degrees = rand(mindegree:maxdegree, nequations)
    vars = map(i -> Symbol("x_$(i)"), 1:nvars)
    randomsystem(T, degrees, vars; rng=rng)
end
function randomsystem(nequations::Int, nvars::Int; kwargs...)
    randomsystem(Complex128, nequations, nvars; kwargs...)
end
function randomsystem(::Type{T}, degrees::Vector{Int}, vars::Vector{Symbol}; rng=Base.Random.GLOBAL_RNG) where {T<:Number}
    map(degrees) do degree
        exponents = create_exponents(degree, length(vars))
        coeffs = rand(rng, T, length(exponents))
        FP.Polynomial(hcat(exponents...), coeffs, vars)
    end
end
function randomsystem(degrees::Vector{Int}, vars::Vector{Symbol}; kwargs...)
    randomsystem(Complex128, degrees, vars; kwargs...)
end

function exponents_helper(curr_sum::Int, target_sum::Int, remaining_elements::Int)::Vector{Vector{Int}}
    if remaining_elements == 0
        return [[]]
    end
    if curr_sum == target_sum
        return [zeros(Int, remaining_elements)]
    end
    if remaining_elements == 1
        return map(x-> [x], 0:(target_sum - curr_sum))
    end

    results = []
    for x=0:(target_sum-curr_sum)
        remaining_results = exponents_helper(curr_sum + x, target_sum, remaining_elements - 1)
        append!(results, map(xs -> [x; xs], remaining_results))
    end
    results
end
create_exponents(total_degree, nvars) = exponents_helper(0, total_degree, nvars)
