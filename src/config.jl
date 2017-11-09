export PolynomialHomotopyConfig, config

"""
    config(H::AbstractHomotopy, [x])

Create a homotopy corresponding to the homotopy type.
For a `AbstractPolynomialHomotopy` this is just `PolynomialHomotopyConfig`.
"""
config(H::AbstractPolynomialHomotopy) = PolynomialHomotopyConfig(H)
config(H::AbstractPolynomialHomotopy, x) = PolynomialHomotopyConfig(H, x)

"""
    PolynomialHomotopyConfig(H::AbstractPolynomialHomotopy{T}, [x::AbstractVector{S}])

A data structure with which `H` and it's derivatives can be evaluated efficiently.
Note that `x` is only used to determine the
output type of `H(x)`.

    PolynomialHomotopyConfig(H::AbstractPolynomialHomotopy{T}, [S])

Instead of a vector `x` a type can also be given directly.
"""
mutable struct PolynomialHomotopyConfig{T} <: AbstractHomotopyConfig{T}
    start::FP.JacobianConfig{T}
    target::FP.JacobianConfig{T}

    result_start::FP.JacobianDiffResult{T, Vector{T}, Matrix{T}}
    result_target::FP.JacobianDiffResult{T, Vector{T}, Matrix{T}}
end

function PolynomialHomotopyConfig(H::AbstractPolynomialHomotopy{T}) where T
    start = FP.JacobianConfig(H.start)
    target = FP.JacobianConfig(H.target)
    result_start = FP.JacobianDiffResult(start)
    result_target = FP.JacobianDiffResult(target)
    PolynomialHomotopyConfig(start, target, result_start, result_target)
end

function PolynomialHomotopyConfig(H::AbstractPolynomialHomotopy, ::AbstractVector{T}) where T
    start = FP.JacobianConfig(H.start, T)
    target = FP.JacobianConfig(H.target, T)
    result_start = FP.JacobianDiffResult(start)
    result_target = FP.JacobianDiffResult(target)
    PolynomialHomotopyConfig(start, target, result_start, result_target)
end

function evaluate_start_target!(
    cfg::PolynomialHomotopyConfig,
    H::AbstractPolynomialHomotopy,
    x::Vector,
    precomputed=false)

    FP.evaluate!(cfg.result_start.value, H.start, x, cfg.start, precomputed)
    FP.evaluate!(cfg.result_target.value, H.target, x, cfg.target, precomputed)
end

function jacobian_start_target!(
    cfg::PolynomialHomotopyConfig,
    H::AbstractPolynomialHomotopy,
    x::Vector,
    precomputed=false)

    FP.jacobian!(cfg.result_start.jacobian, H.start, x, cfg.start, precomputed)
    FP.jacobian!(cfg.result_target.jacobian, H.target, x, cfg.target, precomputed)
end

function evaluate_and_jacobian_start_target!(
    cfg::PolynomialHomotopyConfig,
    H::AbstractPolynomialHomotopy,
    x::Vector)
    FP.jacobian!(cfg.result_start, H.start, x, cfg.start)
    FP.jacobian!(cfg.result_target, H.target, x, cfg.target)
end

value_start(cfg::PolynomialHomotopyConfig) = cfg.result_start.value
value_target(cfg::PolynomialHomotopyConfig) = cfg.result_target.value

jacobian_start(cfg::PolynomialHomotopyConfig) = cfg.result_start.jacobian
jacobian_target(cfg::PolynomialHomotopyConfig) = cfg.result_target.jacobian
