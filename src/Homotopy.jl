__precompile__()

module Homotopy
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

    include("straightline.jl")
    include("geodesic_on_the_sphere.jl")

    include("polynomial.jl")

    include("totaldegree.jl")
    include("randomhomotopy.jl")

    include("condition.jl")
end
