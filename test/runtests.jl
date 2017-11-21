using Homotopies
using Base.Test
import FixedPolynomials
const FP = FixedPolynomials
import DynamicPolynomials: @polyvar

include("straightline_test.jl")
include("geodesic_test.jl")
include("misc_test.jl")
include("condition_test.jl")
