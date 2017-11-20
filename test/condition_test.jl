@testset "Condition StraightLineHomotopy" begin
    @polyvar x
    f = [x^2-1]
    g = [x^2+1]

    z=[1.0]
    t=1.0

    H = StraightLineHomotopy(f,g)
    @test κ(H,z,t) ≈ 1/sqrt(2)
    @test κ_norm(H,z,t) ≈ 1.0
    @test μ_norm(H,z,t) ≈ 1.0
end


@testset "Condition GeodesicOnTheSphere" begin
    @polyvar x
    f = [x^2-1]
    g = [x^2+1]

    z=[1.0, 1.0]
    t=1.0

    H = GeodesicOnTheSphere(f,g)
    @test κ(H,z,t) ≈ 1/sqrt(2)
    @test κ_norm(H,z,t) ≈ 1.0
    @test μ_norm(H,z,t) ≈ 1.0
end
