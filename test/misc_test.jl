@testset "Total degree" begin
    @polyvar x y z

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

    f = x^2+y^2+y*z
    g = x*z + y * z

    H, s = totaldegree(GeodesicOnTheSphere, [f, g])
    @test nvariables(H) == 3

end


@testset "randomsystem" begin
    f = randomsystem(Complex128, 1, 4, density=1.0)[]
    d = FP.degree(f)
    @test length(f.coefficients) == binomial(d + 4, 4)

    f = randomsystem(Complex128, 1, 4, density=0.5)[]
    d = FP.degree(f)
    @test length(f.coefficients) == ceil(Int, binomial(d + 4, 4) * 0.5)

    @test_throws AssertionError randomsystem(Complex128, 1, 4, density=0.0)
    @test_throws AssertionError randomsystem(Complex128, 1, 4, density=-12.0)
    @test_throws AssertionError randomsystem(Complex128, 1, 4, density=213.0)
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

@testset "gammatrick" begin
    @polyvar x y

    H = StraightLineHomotopy{Complex128}([x^3*y-2, x^2+y+1], [x^2+y+1, x^3*y-2])
    G = GeodesicOnTheSphere{Complex128}([x^3*y-2, x^2+y+1], [x^2+y+1, x^3*y-2])
    H2 = deepcopy(H)

    x = rand(2)
    t = 0.535

    @test H(x, t) == H2(x, t)

    gammatrick!(H)
    @test H(x, t) != H2(x, t)

    H = deepcopy(H2)
    H1 = deepcopy(H2)
    gammatrick!(H, 21312)
    gammatrick!(H1, 21312)

    @test H(x, t) == H1(x, t)

    H = deepcopy(H2)
    gammatrick!(H, 0.5)
    @test H(x, 1.0) ≈ 0.5 * H2(x, 1.0)


    G2 = deepcopy(G)
    x = rand(3)
    gammatrick!(G)
    @test G(x, t) != G2(x, t)
end
