using BenchmarkTools



exps = create_exponents(6, 10)

@btime rand($exps, 232)


binomial(4+3, 3)

using Homotopies





import DynamicPolynomials: @polyvar
import FixedPolynomials
const FP = FixedPolynomials

@polyvar x y
f = x^3*y-x*y
g = x+3y*x^2*y+y^3

H = homogenize(GeodesicOnTheSphere{Complex128}([f, g], [2*g, f]))
H2 = homogenize(StraightLineHomotopy{Complex128}([f, g], [2*g, f]))

JH = jacobian!(H)

u = zeros(Complex128, 2)
w = rand(Complex128, 3)
t = rand()
@benchmark evaluate!($u, $H, $w, $t)
@benchmark evaluate!($u, $H2, $w, $t)

v = zeros(Complex128, 2, 2)

@benchmark JH($v, $w, $t)

x = rand(Complex128, 2)

JH = jacobian!(H)

u = zeros(Complex128, 2, 2)

@benchmark JH($u, $x, $one(Complex128))


J_start = [FP.differentiate(f, i) for f in H.start, i=1:FP.nvariables.(H.start[1])]
p = J_start[1, 1]
y = rand(2)

eval2(f) = FP.evaluate(f, x)

@benchmark map!($eval2, $u, $J_start)

@benchmark FP.evaluate($p, $y)
