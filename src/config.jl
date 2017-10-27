export PolynomialConfig

"""
    PolynomialConfig(H::AbstractPolynomialHomotopy{T}, [x::AbstractVector{S}])

A data structure with which `H` and it's derivatives can be evaluated efficiently.
Note that `x` is only used to determine the
output type of `H(x)`.

    PolynomialConfig(H::AbstractPolynomialHomotopy{T}, [S])

Instead of a vector `x` a type can also be given directly.
"""
mutable struct PolynomialConfig{T}
    start::FP.JacobianConfig{T}
    target::FP.JacobianConfig{T}

    result_start::FP.JacobianDiffResult{T}
    result_target::FP.JacobianDiffResult{T}
end

function PolynomialConfig(H::AbstractPolynomialHomotopy{T}) where T
    start = FP.JacobianConfig(H.start)
    target = FP.JacobianConfig(H.target)
    result_start = FP.JacobianDiffResult(start)
    result_target = FP.JacobianDiffResult(target)
    PolynomialConfig(start, target, result_start, result_target)
end

function PolynomialConfig(H::AbstractPolynomialHomotopy, ::AbstractVector{T}) where T
    start = FP.JacobianConfig(H.start, T)
    target = FP.JacobianConfig(H.target, T)
    result_start = FP.JacobianDiffResult(start)
    result_target = FP.JacobianDiffResult(target)
    PolynomialConfig(start, target, result_start, result_target)
end

function evaluate_start_target!(
    cfg::PolynomialConfig,
    H::AbstractPolynomialHomotopy,
    x::Vector)

    FP.evaluate!(cfg.result_start.value, H.start, x, cfg.start)
    FP.evaluate!(cfg.result_target.value, H.target, x, cfg.target)
end

function jacobian_start_target!(
    cfg::PolynomialConfig,
    H::AbstractPolynomialHomotopy,
    x::Vector)

    FP.jacobian!(cfg.result_start.jacobian, H.start, x, cfg.start)
    FP.jacobian!(cfg.result_target.jacobian, H.target, x, cfg.target)
end

function evaluate_and_jacobian_start_target!(
    cfg::PolynomialConfig,
    H::AbstractPolynomialHomotopy,
    x::Vector
    )
    FP.jacobian!(cfg.result_start, H.start, x, cfg.start)
    FP.jacobian!(cfg.result_target, H.target, x, cfg.target)
end

value_start(cfg::PolynomialConfig) = cfg.result_start.value
value_target(cfg::PolynomialConfig) = cfg.result_target.value

jacobian_start(cfg::PolynomialConfig) = cfg.result_start.jacobian
jacobian_target(cfg::PolynomialConfig) = cfg.result_target.jacobian
