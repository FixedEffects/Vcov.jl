@testset "GroupedArray" begin
    N = 10
    a1 = collect(1:N)
    g1 = group(a1)
    @test g1 == GroupedArray(collect(UInt32, 1:N), N)
    @test size(g1) == (N,)
    @test length(g1) == N
    @test g1[1] == UInt32(1)
    @test g1[1:2] == [UInt32(1), UInt32(2)]
    @test g1[g1.<=2] == g1[[1,2]] == g1[1:2]

    a2 = [1,2]
    @test_throws DimensionMismatch group(a1, a2)

    a = rand(N)
    g = group(a)
    @test factorize!(g) == g
end
