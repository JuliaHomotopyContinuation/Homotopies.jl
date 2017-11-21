# Introduction

`Homotopies.jl` is a package for constructing (polynomial) [homotopies](https://en.wikipedia.org/wiki/Homotopy) ``H(x,t)``.

Each implemented homotopy has the same [Interface](@ref) so that you can switch easily between
different homotopy types.
Based on this interface there are also some convenient [higher level constructs](@ref higherlevelconstructs) provided, e.g. the
construction of a total degree system and its start solutions.


## Example
```julia
using Homotopies
# we use an MultivariatePolynomials implementation to construct the homotopy.
import DynamicPolynomials: @polyvar

@polyvar x y z

H = StraightLineHomotopy([x + y^3, x^2*y-2y], [x^3+2, y^3+2])
# H is now StraightLineHomotopy{Int64},
# but let's assume our algorithm uses Complex128, to avoid unnecessary conversions
# it would be better to make
H = StraightLineHomotopy{Complex128}([x + y^3, x^2*y-2y], [x^3+2, y^3+2])

# we can now evaluate H
evaluate(H, rand(Complex128, 2), 0.42)
# or alternatively
H(rand(Complex128, 2), 0.42)
```


## Homotopies

The following homotopies are implemented
### Polynomial homotopies
These are subtypes of `AbstractPolynomialHomotopy`
```@docs
StraightLineHomotopy
GeodesicOnTheSphere
```
