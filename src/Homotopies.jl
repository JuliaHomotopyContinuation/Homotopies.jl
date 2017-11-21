__precompile__()

module Homotopies
    import MultivariatePolynomials
    const MP = MultivariatePolynomials
    import FixedPolynomials
    const FP = FixedPolynomials
    import Base: length, ==

    abstract type AbstractHomotopy{T<:Number} end
    abstract type AbstractPolynomialHomotopy{T<:Number} <: AbstractHomotopy{T} end
    abstract type AbstractHomotopyConfig{T} end
    export AbstractHomotopy, AbstractPolynomialHomotopy, AbstractHomotopyConfig

    include("interface.jl")

    include("config.jl")
    include("diffresult.jl")

    include("homotopies/straightline.jl")
    include("homotopies/geodesic_on_the_sphere.jl")

    include("constructors.jl")

    include("polynomial.jl")

    include("totaldegree.jl")
    include("randomhomotopy.jl")

    include("condition.jl")

    include("bigfloat_utilities.jl")
end
