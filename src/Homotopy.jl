__precompile__()

module Homotopy
    import MultivariatePolynomials
    const MP = MultivariatePolynomials
    import FixedPolynomials
    const FP = FixedPolynomials
    import Base: length, ==

    abstract type AbstractHomotopy{T<:Number} end
    abstract type AbstractPolynomialHomotopy{T<:Number} <: AbstractHomotopy{T} end
    export AbstractHomotopy, AbstractPolynomialHomotopy

    include("interface.jl")

    include("straightline.jl")
    include("geodesic_on_the_sphere.jl")
    include("gammatrick.jl")

    include("totaldegree.jl")
    include("randomhomotopy.jl")

    include("condition_number/condition_number.jl")
end
