using Homotopy
using Base.Test
import FixedPolynomials
const FP = FixedPolynomials
import DynamicPolynomials: @polyvar

@testset "StraightLineHomotopy" begin
    @polyvar x y z
    f = x + y + 2z
    g = x + z
    h = y + z


    @test StraightLineHomotopy(f, g) isa StraightLineHomotopy{Int64}
    @test StraightLineHomotopy([f, f], [g, h]) isa StraightLineHomotopy{Int64}
    @test StraightLineHomotopy{Float64}([f, g], [f, h]) isa StraightLineHomotopy{Float64}
    @test StraightLineHomotopy{Float64}(f, g) isa StraightLineHomotopy{Float64}
    @test StraightLineHomotopy([FP.Polynomial(f)], [FP.Polynomial{Complex128}(f)]) isa StraightLineHomotopy{Complex128}
    @test_throws AssertionError StraightLineHomotopy{Complex128}([f], [g])
    @test_throws AssertionError StraightLineHomotopy{Complex128}([f], [f, f])


    H = StraightLineHomotopy{Float64}([f, g], [f, h])
    K = StraightLineHomotopy([f, g], [f, h])
    @test promote_type(typeof(H), typeof(K)) == StraightLineHomotopy{Float64}
    @test convert(StraightLineHomotopy{Float64}, K) isa StraightLineHomotopy{Float64}


    w = [1.0, 2.0, 2.0]
    @test evaluate(H, w, 1.0) == [7, 3]
    @test H(w, 1.0) == [7, 3]
    u = zeros(2)
    evaluate!(u, H, w, 1.0)
    @test u == [7, 3]


    JH = jacobian(H)
    @test JH(rand(3), 0.0) == [1 1 2; 0 1 1]
    @test JH(rand(3), 1.0) == [1 1 2; 1 0 1]
    u = zeros(2, 3)
    JH! = jacobian!(H)
    @test JH!(u, rand(3), 0.0) == JH(rand(3), 0.0)
    @test u == JH(rand(3), 0.0)


    HDT = dt(H)
    @test HDT([1.0, 2, 2], 1.0) == [0, -1]
    u = zeros(2)
    HDT! = dt!(H)
    HDT!(u, [1, 2, 2.0], 1.0)
    @test u == [0, -1]
    @test string(H) == "Homotopy.StraightLineHomotopy{Float64}((1-t)⋅[x+y+2.0z, y+z] + t⋅[x+y+2.0z, x+z])\n"

    N=weylnorm(H)
    @test N(0.0)==sqrt(8)
    @test N(0.5)==sqrt(7.5)

    H = StraightLineHomotopy(x^2+y+z, z^2+2+x+y^2)
    @test ishomogenous(H) == false
    @test ishomogenized(H) == false
    HH = homogenize(H)
    @test ishomogenous(HH)
    @test ishomogenized(HH)
    @test H != HH
    @test HH == homogenize(HH)
    @test isequal(HH, homogenize(HH))
    @test dehomogenize(HH) == H

    @test length(H) == 1
    @test nvariables(H) == 3
end

@testset "GammaTrickHomotopy" begin
    @polyvar x y z
    f = x + y + 2z
    g = x + z
    h = y + z

    γ = rand(Complex128)
    @test GammaTrickHomotopy(f, g) isa GammaTrickHomotopy{Complex128}
    @test GammaTrickHomotopy([f, f], [g, h]) isa GammaTrickHomotopy{Complex128}
    @test GammaTrickHomotopy([f, f], [g, h], 2323) isa GammaTrickHomotopy{Complex128}
    @test GammaTrickHomotopy([f, f], [g, h], γ) isa GammaTrickHomotopy{Complex128}
    @test GammaTrickHomotopy{Complex64}([f, g], [f, h]) isa GammaTrickHomotopy{Complex64}
    @test GammaTrickHomotopy{Complex64}([f, g], [f, h], 232) isa GammaTrickHomotopy{Complex64}
    @test GammaTrickHomotopy{Complex64}([f, g], [f, h], γ) isa GammaTrickHomotopy{Complex64}
    @test GammaTrickHomotopy{Complex64}(f, g) isa GammaTrickHomotopy{Complex64}
    @test GammaTrickHomotopy{Complex64}(f, g, 123) isa GammaTrickHomotopy{Complex64}
    @test GammaTrickHomotopy{Complex64}(f, g, γ) isa GammaTrickHomotopy{Complex64}
    @test GammaTrickHomotopy{Complex64}([FP.Polynomial(f)], [FP.Polynomial{Complex128}(f)]) isa GammaTrickHomotopy{Complex64}
    @test GammaTrickHomotopy{Complex64}([FP.Polynomial(f)], [FP.Polynomial{Complex128}(f)], 232) isa GammaTrickHomotopy{Complex64}
    @test GammaTrickHomotopy{Complex64}([FP.Polynomial(f)], [FP.Polynomial{Complex128}(f)], γ) isa GammaTrickHomotopy{Complex64}
    @test GammaTrickHomotopy([FP.Polynomial(f)], [FP.Polynomial{Complex128}(f)]) isa GammaTrickHomotopy{Complex128}
    @test GammaTrickHomotopy([FP.Polynomial(f)], [FP.Polynomial{Complex128}(f)], 12312) isa GammaTrickHomotopy{Complex128}
    @test GammaTrickHomotopy([FP.Polynomial(f)], [FP.Polynomial{Complex128}(f)], γ) isa GammaTrickHomotopy{Complex128}

    @test_throws AssertionError GammaTrickHomotopy{Complex128}([f], [g])
    @test_throws AssertionError GammaTrickHomotopy{Complex128}([f], [f, f])

    H = GammaTrickHomotopy([f, g], [f, h])
    K = GammaTrickHomotopy{Complex32}([f, g], [f, h])
    @test promote_type(typeof(H), typeof(K)) == GammaTrickHomotopy{Complex128}
    @test convert(GammaTrickHomotopy{Complex128}, K) isa GammaTrickHomotopy{Complex128}


    @test GammaTrickHomotopy([f, f], [g, h], 2323).γ == GammaTrickHomotopy([f, f], [g, h], 2323).γ


    w = complex.([1.0, 2.0, 2.0])
    @test evaluate(H, w, 1.0) == H.γ * [7, 3]
    @test H(w, 1.0) == H.γ * [7, 3]
    u = zeros(Complex128, 2)
    evaluate!(u, H, w, 1.0)
    @test u == H.γ * [7, 3]

    JH = jacobian(H)
    @test JH(rand(Complex128, 3), 0.0) == [1 1 2; 0 1 1]
    @test JH(rand(Complex128, 3), 1.0) == H.γ * [1 1 2; 1 0 1]
    u = zeros(Complex128, 2, 3)
    JH! = jacobian!(H)
    @test JH!(u, rand(Complex128, 3), 0.0) == JH(rand(Complex128, 3), 0.0)
    @test u == JH(rand(Complex128, 3), 0.0)


    HDT = dt(H)
    @test HDT([1.0+0im, 2, 2], 1.0) == H.γ * [7, 3] - [7, 4]
    u = zeros(Complex128, 2)
    HDT! = dt!(H)
    HDT!(u, [1+0im, 2, 2.0], 1.0)
    @test u == H.γ * [7, 3] - [7, 4]
    @test string(H) == "Homotopy.GammaTrickHomotopy{Complex{Float64}}((1-t)⋅[x+y+2.0z, y+z] + t⋅$(H.γ)⋅[x+y+2.0z, x+z])\n"

    N=weylnorm(H)
    @test N(0.0)==sqrt(8)


    H = GammaTrickHomotopy(x^2+y+z, z^2+2+x+y^2)
    @test ishomogenous(H) == false
    @test ishomogenized(H) == false
    HH = homogenize(H)
    @test ishomogenous(HH)
    @test ishomogenized(HH)
    @test H != HH
    @test HH == homogenize(HH)
    @test isequal(HH, homogenize(HH))
    @test dehomogenize(HH) == H

    @test length(H) == 1
    @test nvariables(H) == 3
end


@testset "GeodesicOnTheSphere" begin
    @polyvar x y z
    f = x + y + 3z
    g = x + 2z
    h = y + 2z


    @test GeodesicOnTheSphere(f, g) isa GeodesicOnTheSphere{Float64}
    @test GeodesicOnTheSphere([f, f], [g, h]) isa GeodesicOnTheSphere{Float64}
    @test GeodesicOnTheSphere{Float64}([f, g], [f, h]) isa GeodesicOnTheSphere{Float64}
    @test GeodesicOnTheSphere{Float64}(f, g) isa GeodesicOnTheSphere{Float64}
    @test GeodesicOnTheSphere([FP.Polynomial(f)], [FP.Polynomial{Complex128}(f)]) isa GeodesicOnTheSphere{Complex128}
    @test_throws AssertionError GeodesicOnTheSphere{Complex128}([f], [g])
    @test_throws AssertionError GeodesicOnTheSphere{Complex128}([f], [f, f])


    H = GeodesicOnTheSphere{Float64}([f, g], [f, h])
    K = GeodesicOnTheSphere([f, g], [f, h])
    @test promote_type(typeof(H), typeof(K)) == GeodesicOnTheSphere{Float64}
    @test convert(GeodesicOnTheSphere{Float64}, K) isa GeodesicOnTheSphere{Float64}


    w = [1.0, 2.0, 2.0]
    @test evaluate(H, w, 1.0) == [9/4, 5/4]
    @test H(w, 1.0) == [9/4, 5/4]
    u = zeros(2)
    evaluate!(u, H, w, 1.0)
    @test u == [9/4, 5/4]


    JH = jacobian(H)
    @test JH(rand(3), 0.0) == [1 1 3; 0 1 2]./4
    @test JH(rand(3), 1.0) == [1 1 3; 1 0 2]./4
    u = zeros(2, 3)
    JH! = jacobian!(H)
    @test JH!(u, rand(3), 0.0) == JH(rand(3), 0.0)
    @test u == JH(rand(3), 0.0)


    HDT = dt(H)
    @test HDT([1, 2, 2.0], 0.0) == [0, 0]
    u = zeros(2)
    HDT! = dt!(H)
    HDT!(u, [1, 2, 2.0], 0.0)
    @test u == [0, 0]
    @test string(H) == "Homotopy.GeodesicOnTheSphere{Float64} The homotopy is given by the spherical geodesic from start/|start| (t=1) to target/|target| (t=0). \n"

    N=weylnorm(H)
    @test N(0.0)==1.0
    @test N(0.5)==1.0

    H = GeodesicOnTheSphere(x^2+y+z, z^2+2+x+y^2)
    HH = homogenize(H)
    @test ishomogenous(H) == false
    @test ishomogenized(H) == false
    @test ishomogenous(HH)
    @test ishomogenized(HH)
    @test H != HH


    @test length(H) == 1
    @test nvariables(H) == 3
end


@testset "Total degree" begin
    @polyvar x y

    H, solutions = totaldegree(StraightLineHomotopy{Complex128}, [x^2+y+1, x^3*y-2])
    @test length(solutions) == 8
    @test eltype(solutions) == Vector{Complex128}
    for sol in solutions
        @test norm(H(sol, 1.0)) ≈ 0 atol=1e-14
    end
    H, solutions = totaldegree(StraightLineHomotopy{Complex128}, [x^2+y+1, x^3*y-2], unitroots=true)
    for sol in solutions
        @test norm(sol) ≈ norm(ones(2)) atol=1e-14
    end

    H, solutions = totaldegree(GammaTrickHomotopy{Complex128}, [x^2+y+1, x^3*y-2])
    for sol in solutions
        @test norm(H(sol, 1.0)) ≈ 0 atol=1e-14
    end

    H, sols = totaldegree(GammaTrickHomotopy, [x^2+y+1, x^3*y-2])
    @test H isa GammaTrickHomotopy{Complex128}
end


@testset "randomhomotopy" begin
    H, solutions = randomhomotopy(StraightLineHomotopy{Complex64}, 3)
    @test H isa StraightLineHomotopy{Complex64}

    H, solutions = randomhomotopy(StraightLineHomotopy, 3)
    @test H isa StraightLineHomotopy{Complex128}
    @test length(H) == 3
    @test nvariables(H) == 3
    @test length(solutions) == prod(FP.degree.(H.start))

    for sol in solutions
        @test norm(H(sol, 1.0)) ≈ 0 atol=1e-14
    end
end
