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

    w = [1.0, 1.0, 2.0, 2.0]
    @test evaluate(H, w, 1.0) == [9/4, 5/4]
    @test H(w, 1.0) == [9/4, 5/4]
    u = zeros(2)
    evaluate!(u, H, w, 1.0)
    @test u == [9/4, 5/4]

    cfg = PolynomialConfig(H)
    @test evaluate(H, w, 1.0, cfg) == [9/4, 5/4]
    u = zeros(2)
    evaluate!(u, H, w, 1.0, cfg)
    @test u == [9/4, 5/4]

    @test jacobian(H, rand(4), 0.0, cfg) == [0 1 1 3; 0 0 1 2]./4
    @test jacobian(H, rand(4), 1.0, cfg) == [0 1 1 3; 0 1 0 2]./4
    u = zeros(2, 4)
    @test jacobian!(u, H, rand(4), 0.0, cfg) == [0 1 1 3; 0 0 1 2]./4
    @test u == [0 1 1 3; 0 0 1 2]./4


    # TODO: This test is bonkers atm
    @test dt(H, [1.0, 1.0, 1, 2, 2.0], 0.0, cfg) != [0, 0]
    u = zeros(2)
    dt!(u, H, [1.0, 1.0, 1, 2, 2.0], 0.0, cfg)
    @test u == dt(H, [1.0, 1.0, 1, 2, 2.0], 0.0, cfg)
    @test string(H) == "Homotopy.GeodesicOnTheSphere{Float64}((1-t)⋅[0.25x+0.25y+0.75z, 0.25y+0.5z] + t⋅[0.25x+0.25y+0.75z, 0.25x+0.5z]) with angle $(H.α)\n"

    N=weylnorm(H)
    @test N(0.0)==1.0
    @test N(0.5)==1.0

    H = GeodesicOnTheSphere(x^2+y+z, z^2+2+x+y^2)
    #HH = homogenize(H)
    @test ishomogenous(H) == true
    @test ishomogenized(H) == true
    # @test ishomogenous(HH)
    # @test ishomogenized(HH)
    # @test H != HH


    @test length(H) == 1
    @test nvariables(H) == 4
end
