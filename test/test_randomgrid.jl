
function test_randomgrids()
    println("Random grids:")
    @testset begin
        d = UnitDisk()
        g = randomgrid(d, 10)
        @test length(g) == 10
        @test length(eltype(g)) == dimension(d)
        @test reduce(&, [x ∈ d for x in g])

        g2 = randomgrid(UnitDisk{BigFloat}(), 10)
        @test eltype(g2[1]) == BigFloat

        g3 = randomgrid(0.0..1.0, 10)
        @test length(g3) == 10
        # 1D is a special case where we don't use vectors of length 1
        @test eltype(g3) == Float64

        box1 = UnitInterval()^2
        p1 = randompoint(box1)
        @test length(p1) == dimension(box1)
        @test p1 ∈ box1
        box2 = 0.0..1.0
        p2 = randompoint(box2)
        @test typeof(p2) == Float64
        @test p2 ∈ box2

        g3 = randompoint(UnionDomain(0.0..1.5,1.5..2.0))
        @test typeof(g3) == Float64
        @test g3 ∈ 0.0..2.0

        g3 = randompoint(IntersectionDomain(0.0..1.5,1.0..2.0))
        @test typeof(g3) == Float64
        @test g3 ∈ 1.0..1.5

        g3 = randompoint(DifferenceDomain(0.0..1.5,1.0..2.0))
        @test typeof(g3) == Float64
        @test g3 ∈ 0.0..1.0
    end
end

test_randomgrids()
