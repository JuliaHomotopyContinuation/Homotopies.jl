export evaluate, evaluate!, jacobian, jacobian!, dt, dt!,
    nvariables, homogenize, dehomogenize, ishomogenized, ishomogenous

"""
    evaluate(H::AbstractHomotopy, x, t)

Evaluate the homotopy `H` at `x` to time `t`, i.e. ``H(x,t)``
"""
function evaluate end

"""
    evaluate!(u::Vector, H::AbstractHomotopy, x, t)

Evaluate the homotopy `H` at `x` to time `t`, i.e. `H(x,t)`, and store the result in `u`.
Use this instead of [`evaluate`](@ref) to avoid allocations.
"""
function evaluate! end

"""
    jacobian(H::AbstractHomotopy)

Compute an evaluation function `(x, t) -> J_H(x,t)` of the jacobian ``J_H`` of the homotopy ``H``.
The jacobian is constructed w.r.t. `x`, i.e. it doesn't contain the partial derivatives
w.r.t. `t`.
"""
function jacobian end

"""
    jacobian!(H::AbstractHomotopy)

Compute an inplace evaluation function `(u, x, t) -> u := J_H(x,t)` of the jacobian
``J_H`` of the homotopy ``H``.
Use this instead of [`jacobian`](@ref) to avoid allocations.
"""
function jacobian! end

"""
    dt(H::AbstractHomotopy)

Compute an evaluation function `(x, t) -> ∂H∂t(x,t)` of the partial derivative
``\frac{∂H}{∂t}`` of the homotopy ``H``.
"""
function dt end

"""
    dt!(H::AbstractHomotopy)

Compute an inplace evaluation function `(u, x, t) -> u := ∂H∂t(x,t)` of the partial derivative
``\frac{∂H}{∂t}`` of the homotopy ``H``. Use this instead of [`dt`](@ref) to avoid allocations.
"""
function dt! end

"""
    nvariables(H::AbstractHomotopy)

The number of variables which `H` expects as input, i.e. to evaluate `H(x,t)` `x` has to be a
vector of length `nvariables(H)`.
"""
function nvariables end

"""
    homogenize(H::AbstractHomotopy)

Homogenize the homotopy `H`. This adds an additional variable.
If `H` is already homogenized, this is the identity.
"""
function homogenize end

"""
    dehomogenize(H::AbstractHomotopy)

Dehomogenize the homotopy `H`. This removes the first variable.
If `H` is not homogenized, this is the identity.
"""
function dehomogenize end

"""
    ishomogenized(H::AbstractHomotopy)

Check whether the homotopy `H` was homogenized.
"""
function ishomogenized end

"""
    ishomogenous(H::AbstractHomotopy)

Check whether the homotopy `H` is homogenous. This does not imply that `H` was homogenized.
"""
function ishomogenous end
