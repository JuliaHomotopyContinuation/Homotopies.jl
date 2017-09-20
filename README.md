# Homotopy

A package for constructing [homotopies](https://en.wikipedia.org/wiki/Homotopy) `H(x,t)`.
Currently the following homotopies are implemented:

**StraightLineHomotopy**

A `StraightLineHomotopy` has the form
```julia
(1 - t) * F + t * G
```
where `F` and `G` are systems of polynomials.


**GammaTrickHomotopy**


A `GammaTrickHomotopy` has the form
```julia
(1 - t) * F + t * γ * G`
```
where `F` and `G` are systems of polynomials and `γ` is random complex number.

The systems of polynomials are represented as `Vector`s of [`FixedPolynomials.Polynomial`](https://github.com/saschatimme/FixedPolynomials.jl)s.
