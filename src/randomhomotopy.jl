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
    randomsystem([T=Complex128,] nequations::Int, nvars::Int; mindegree=0, maxdegree=5, rng=Base.Random.GLOBAL_RNG, density=rand() * 0.8 + 0.1)

Creates a random polynomial system of `nequations` equations with `nvars` variables (named ``x_1``, ...``x_{nvars}``).
Each polynomial has a total degree uniformly drawn from ``[mindegree, maxdegree]``.
The coefficients are drawn independently from the given `rng`.
With `density` you can control how many coefficients are non-zero. A value of `1.0` creates
a dense polynomial (i.e. every coefficient is non-zero). A value of `0.5` creates a polynomial
where only half of all monomials are non zero.

    randomsystem([T=Complex128,] degrees::Vector{Int}, variables::Vector{Symbol}; rng=N(0,1))

Create a random polynomial system with the given `degrees` and `variables`.
"""
function randomsystem(::Type{T}, nequations::Int, nvars::Int; mindegree=2, maxdegree=5, rng=Base.Random.GLOBAL_RNG, density=(rand() * 0.8 + 0.1)) where {T<:Number}
    @assert (0 < density ≤ 1) "Expexted a density between 0.0 and 1.0"
    degrees = rand(mindegree:maxdegree, nequations)
    vars = map(i -> Symbol("x_$(i)"), 1:nvars)
    randomsystem(T, degrees, vars, rng, density)
end
function randomsystem(nequations::Int, nvars::Int; kwargs...)
    randomsystem(Complex128, nequations, nvars; kwargs...)
end
function randomsystem(::Type{T}, degrees::Vector{Int}, vars::Vector{Symbol}, rng, density) where {T<:Number}
    map(degrees) do degree
        exponents = create_exponents(degree, length(vars))
        nexponents = length(exponents)
        indices = rand(1:nexponents, ceil(Int, nexponents * density))
        exp_matrix = Matrix{Int}(length(vars), length(indices))
        for (j, index) in enumerate(indices)
            exp_matrix[:,j] .= exponents[index]
        end
        coeffs = randomcoeffs(rng, T, length(indices))
        FP.Polynomial(exp_matrix, coeffs, vars)
    end
end

randomcoeffs(rng, ::Type{Complex{BigFloat}}, n) = convert.(Complex{BigFloat}, rand(rng, Complex128, n))
randomcoeffs(rng, ::Type{BigFloat}, n) = convert.(BigFloat, rand(rng, Float64, n))
randomcoeffs(rng, ::Type{T}, n) where {T<:Number} = rand(rng, T, n)

function randomsystem(degrees::Vector{Int}, vars::Vector{Symbol}; kwargs...)
    randomsystem(Complex128, degrees, vars; kwargs...)
end

function exponents_helper!(results, partial_result, curr_sum::Int, target_sum::Int, remaining_elements::Int)::Vector{Vector{Int}}
    if remaining_elements == 0
        push!(results, partial_result)
        return results
    end
    if curr_sum == target_sum
        for _ = 1:remaining_elements
            push!(partial_result, 0)
        end
        push!(results, partial_result)
        return results
    end
    if remaining_elements == 1
        for x = 0:(target_sum - curr_sum) - 1
            push!(results, [partial_result; x])
        end
        push!(partial_result, target_sum - curr_sum)
        push!(results, partial_result)
        return results
    end

    for x=0:(target_sum-curr_sum)
        exponents_helper!(results, [partial_result; x], curr_sum + x, target_sum, remaining_elements - 1)
    end
    results
end

create_exponents(total_degree, nvars) = exponents_helper!(Vector{Vector{Int}}(), Int[], 0, total_degree, nvars)
