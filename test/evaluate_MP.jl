@testset "Evaluate MultivariatePolynomials" begin
    @polyvar x
    f = [x^2-1]
    g = [x^2+1]

    z = [1.0]
    u = [1]
    v = [1.0 + im * 0.0]


    @test evaluate(f, z) == [0.0]
    @test evaluate(g, z) == [2.0]

    @test evaluate(f, u) == [0]
    @test evaluate(g, u) == [2]

    @test evaluate(f, v) == [0.0 + im * 0.0]
    @test evaluate(g, v) == [2.0 + im * 0.0]
end
