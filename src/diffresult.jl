export DtDiffResult, value, dt, JacobianDiffResult

"""
    DtDiffResult(cfg::PolynomialHomotopyConfig)

During the computation of the time derivative ``∂H/∂t`` we compute nearly everything
we need for the evaluation of ``H(x)``. `DtDiffResult` allocates memory to hold both values.
This structure also signals `dt!` to store ``H(x)`` and ``∂H/∂t``.

### Example

```julia
cfg = PolynomialHomotopyConfig(H, x)
r = JacobianDiffResult(cfg)
dt!(r, H, x, t, cfg)

value(r) == evaluate(H, x, t, cfg)
dt(r) == dt(H, x, t, cfg)
```

    DtDiffResult(value::AbstractVector, jacobian::AbstractMatrix)

Allocate the memory to hold the value and the derivative ``∂H/∂t`` by yourself.
"""
mutable struct DtDiffResult{T, AV1<:AbstractVector{T}, AV2<:AbstractVector{T}}
    value::AV1
    dt::AV2
end

function DtDiffResult(cfg::PolynomialHomotopyConfig{T}) where T
    DtDiffResult{T, Vector{T}, Vector{T}}(
        similar(cfg.result_start.value),
        similar(cfg.result_start.value))
end

function DtDiffResult(value::AbstractVector{T}, dt::AbstractVector{T}) where T
    DtDiffResult{T, typeof(value), typeof(jacobian)}(value, dt)
end
"""
    value(r::DtDiffResult)

Get the value stored in `r`.
"""
value(r::DtDiffResult) = r.value


"""
    dt(r::DtDiffResult)

Get the time derivative stored in `r`.
"""
dt(r::DtDiffResult) = r.dt


"""
    JacobianDiffResult(cfg::PolynomialHomotopyConfig)

During the computation of the jacobian ``J_H(x)`` we compute nearly everything we need for the evaluation of
``H(x)``. `JacobianDiffResult` allocates memory to hold both values.
This structure also signals `jacobian!` to store ``H(x)`` and ``J_H(x)``.

### Example

```julia
cfg = PolynomialHomotopyConfig(H, x)
r = JacobianDiffResult(cfg)
jacobian!(r, H, x, t, cfg)

value(r) == evaluate(H, x, t, cfg)
jacobian(r) == jacobian(H, x, t, cfg)
```

    JacobianDiffResult(value::AbstractVector, jacobian::AbstractMatrix)

Allocate the memory to hold the value and the jacobian by yourself.
"""
mutable struct JacobianDiffResult{T, AV<:AbstractVector{T}, AM<:AbstractMatrix{T}}
    value::AV
    jacobian::AM
end

function JacobianDiffResult(cfg::PolynomialHomotopyConfig{T}) where T
    JacobianDiffResult{T, Vector{T}, Matrix{T}}(
        similar(cfg.result_start.value),
        similar(cfg.result_start.jacobian))
end

function JacobianDiffResult(value::AbstractVector{T}, jacobian::AbstractMatrix{T}) where T
    JacobianDiffResult{T, typeof(value), typeof(jacobian)}(value, jacobian)
end

"""
    value(r::JacobianDiffResult)

Get the value stored in `r`.
"""
value(r::JacobianDiffResult) = r.value

"""
    jacobian(r::JacobianDiffResult)

Get the jacobian stored in `r`.
"""
jacobian(r::JacobianDiffResult) = r.jacobian
