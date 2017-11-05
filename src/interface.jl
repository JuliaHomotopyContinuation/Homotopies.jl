export evaluate, evaluate!, weylnorm, jacobian, jacobian!, dt, dt!,
    nvariables, homogenize, dehomogenize, ishomogenized, ishomogenous

"""
    evaluate(H::AbstractPolynomialHomotopy, x, t)

Evaluate the homotopy `H` at `x` to time `t`, i.e. ``H(x,t)``.

    evaluate(H::AbstractPolynomialHomotopy, x, t, cfg::PolynomialHomotopyConfig)

Evaluate the homotopy `H` at `x` to time `t` using the precompuated values in `cfg`.
Note that this is significantly faster than `evaluate(H, x, t)`.
"""
function evaluate end

"""
    evaluate!(u::Vector, H::AbstractPolynomialHomotopy, x, t)

Evaluate the homotopy `H` at `x` to time `t`, i.e. ``H(x,t)``, and store the result in `u`.

    evaluate!(u::AbstractVector, H::AbstractPolynomialHomotopy, x, t, cfg::PolynomialHomotopyConfig)

Evaluate the homotopy `H` at `x` to time `t` using the precompuated values in `cfg` and store
the result in `u`.
"""
function evaluate! end


"""
    weylnorm(H::AbstractPolynomialHomotopy)

Creates a function with variable `t` that computes the Weyl norm (or Bombieri norm) of ``H(x,t)``.
See [here](https://en.wikipedia.org/wiki/Bombieri_norm) for details about the Weyl norm.
"""
function weylnorm end

"""
    jacobian(H::AbstractPolynomialHomotopy, x, t, cfg::PolynomialHomotopyConfig)

Compute the jacobian of `H` at `x` and `t` using the precomputed values in `cfg`.
The jacobian is constructed w.r.t. `x`, i.e. it doesn't contain the partial derivatives
w.r.t. `t`.
"""
function jacobian end

"""
    jacobian!(u, H::AbstractPolynomialHomotopy, x, t, cfg::PolynomialHomotopyConfig)

Compute the jacobian of `H` at `x` and `t` using the precomputed values in `cfg` and
store the result in `u`.

    jacobian!(r::JacobianDiffResult, H::AbstractPolynomialHomotopy, x, t, cfg::PolynomialHomotopyConfig)

Compute ``H(x, t)`` and the jacobian of `H` at `x` and `t` at once using the precomputated values in `cfg`
and store thre result in `r`. This is faster than computing both values separetely.

### Example
```julia
cfg = PolynomialHomotopyConfig(H)
r = JacobianDiffResult(cfg)
jacobian!(r, H, x, t, cfg)

value(r) == H(x, t)
jacobian(r) == jacobian(H, x, t, cfg)
```
"""
function jacobian! end

"""
    dt(H::AbstractPolynomialHomotopy, x, t, cfg::PolynomialHomotopyConfig)

Compute the derivative of `H` w.r.t. ``t`` at `x` and `t` using the precomputed values in `cfg`.
"""
function dt end

"""
    dt!(u, H::AbstractPolynomialHomotopy, x, t, cfg::PolynomialHomotopyConfig)

Compute the derivative of `H` w.r.t. ``t`` at `x` and `t` using the precomputed values in `cfg`
and store the result in `u`.

    dt!(r::DtDiffResult, H::AbstractPolynomialHomotopy, x, t, cfg::PolynomialHomotopyConfig)

Compute the derivative of `H` w.r.t. ``t`` at `x` and `t` using the precomputed values in `cfg`
and store the result in `r`. This is faster than computing both values separetely.

### Example
```julia
cfg = PolynomialHomotopyConfig(H)
r = DtDiffResult(cfg)
dt!(r, H, x, t, cfg)

value(r) == H(x, t)
dt(r) == dt(H, x, t, cfg)
```
"""
function dt! end

"""
    nvariables(H::AbstractPolynomialHomotopy)

The number of variables which `H` expects as input, i.e. to evaluate `H(x,t)` `x` has to be a
vector of length `nvariables(H)`.
"""
function nvariables end

"""
    homogenize(H::AbstractPolynomialHomotopy)

Homogenize the homotopy `H`. This adds an additional variable.
If `H` is already homogenized, this is the identity.
"""
function homogenize end

"""
    dehomogenize(H::AbstractPolynomialHomotopy)

Dehomogenize the homotopy `H`. This removes the first variable.
If `H` is not homogenized, this is the identity.
"""
function dehomogenize end

"""
    ishomogenized(H::AbstractPolynomialHomotopy)

Check whether the homotopy `H` was homogenized.
"""
function ishomogenized end

"""
    ishomogenous(H::AbstractPolynomialHomotopy)

Check whether the homotopy `H` is homogenous. This does not imply that `H` was homogenized.
"""
function ishomogenous end
