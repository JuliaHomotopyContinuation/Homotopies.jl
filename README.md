# Homotopy

| **Documentation** | **Build Status** |
|:-----------------:|:----------------:|
| [![][docs-stable-img]][docs-stable-url] | [![Build Status][build-img]][build-url] |
| [![][docs-latest-img]][docs-latest-url] | [![Codecov branch][codecov-img]][codecov-url] |


A package for constructing (polynomial) [homotopies](https://en.wikipedia.org/wiki/Homotopy) `H(x,t)`.
Currently the following homotopies are implemented:

**StraightLineHomotopy**

A `StraightLineHomotopy` has the form
```julia
(1 - t) * F + t * G
```

**GammaTrickHomotopy**

A `GammaTrickHomotopy` has the form
```julia
(1 - t) * F + t * γ * G`
```
where `γ` is random complex number.

## Getting started

You can install the package via
```julia
Pkg.add("https://github.com/JuliaHomotopyContinuation/Homotopy.jl.git")
```

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-url]: https://JuliaHomotopyContinuation.github.io/Homotopy.jl/stable
[docs-latest-url]: https://JuliaHomotopyContinuation.github.io/Homotopy.jl/latest

[build-img]: https://travis-ci.org/JuliaHomotopyContinuation/Homotopy.jl.svg?branch=master
[build-url]: https://travis-ci.org/JuliaHomotopyContinuation/Homotopy.jl
[codecov-img]: https://codecov.io/gh/juliahomotopycontinuation/Homotopy.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/juliahomotopycontinuation/Homotopy.jl
