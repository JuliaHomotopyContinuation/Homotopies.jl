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

    cfg = PolynomialConfig(H)

    @test evaluate(H, w, 1.0, cfg) == [7, 3]
    u = zeros(2)
    evaluate!(u, H, w, 1.0, cfg)
    @test u == [7, 3]



    @test jacobian(H, rand(3), 0.0, cfg) == [1 1 2; 0 1 1]
    @test jacobian(H, rand(3), 1.0, cfg) == [1 1 2; 1 0 1]
    u = zeros(2, 3)
    @test jacobian!(u, H, rand(3), 0.0, cfg) == [1 1 2; 0 1 1]
    @test u == [1 1 2; 0 1 1]


    @test dt(H, [1.0, 2, 2], 1.0, cfg) == [0, -1]
    u = zeros(2)
    dt!(u, H, [1, 2, 2.0], 1.0, cfg)
    @test u == [0, -1]
    @test string(H) == "Homotopy.StraightLineHomotopy{Float64}((1-t)⋅[x+y+2.0z, y+z] + t⋅[x+y+2.0z, x+z])\n"

    N = weylnorm(H)
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
