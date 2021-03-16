
function test_generic_subgrid(grid, s)
    @test covering(grid) ≈ s
    for x in grid
        @test x ∈ covering(grid)
    end
    for x in subindices(grid)
        @test issubindex(x, grid)
    end

    cnt = 0
    for i in eachindex(supergrid(grid))
        if issubindex(i, grid)
            cnt += 1
        end
    end
    @test cnt == length(grid)
end

delimit("Grid functionality")
@testset  "SubGrids" begin
    n = 20
    grid1 = EquispacedGrid(n, -1.0, 1.0)
    subgrid1 = MaskedGrid(grid1, -0.5..0.7)
    test_generic_subgrid(subgrid1, -0.5..0.7)

    subgrid2 = IndexSubGrid(grid1, 4:12)
    subgrid3 = subgrid(grid1, -0.5..0.7)
    test_generic_subgrid(subgrid3, -0.5..0.7)

    @test subgrid1 == subgrid3
    @test mask(subgrid1) == mask(subgrid3)
    @test subindices(subgrid1) == subindices(subgrid3)

    @test covering(subgrid1) ⊆ covering(grid1)
    @test covering(subgrid2) ⊆ covering(grid1)

    G1 = EquispacedGrid(n, -1.0, 1.0)
    G2 = EquispacedGrid(n, -1.0, 1.0)
    ProductG = G1 × G2

    C = UnitDisk()
    refgrid = MaskedGrid(ProductG, C)
    circle_grid = subgrid(ProductG, C)
    @test circle_grid isa MaskedGrid
    @test refgrid ≈ circle_grid
    test_generic_subgrid(circle_grid, C)
    @testset "Generic MaskedGrid" begin

        @test (length(circle_grid)/length(ProductG)-pi*0.25) < 0.01

        G1s = IndexSubGrid(G1,2:4)
        G2s = IndexSubGrid(G2,3:5)
        ProductGs = G1s × G2s
        @test G1s[1] == G1[2]
        @test G2s[1] == G2[3]
        @test ProductGs[1,1] == [G1[2],G2[3]]
    end

    # subgrid of a subgrid
    g1 = EquispacedGrid(10,-1,1)^2
    @test mask(subgrid(subgrid(g1,(0.0..1.0)^2),UnitSimplex{2}()))==[i+j<6 for i in 1:5 , j in 1:5]

    g1 = EquispacedGrid(10,0,1)^2
    for x in g1[BitArray([i+j<12 for i in 1:10 , j in 1:10])]
        @test x ∈ UnitSimplex{2}()
    end

    C = UnitInterval()^2
    productgrid = subgrid(ProductG, C)
    refgrid = MaskedGrid(ProductG, C)

    @test supergrid(productgrid) == ProductG
    @test productgrid isa TensorSubGrid
    test_generic_subgrid(productgrid, C)
    refgrid = MaskedGrid(ProductG, C)
    @test reshape(refgrid,10,10) == productgrid
    @test subindices(refgrid) == subindices(productgrid)


    # Generic tests for the subgrids
    @testset "result" for (grid,subgrid) in ( (grid1,subgrid1), (grid1,subgrid2), (ProductG, circle_grid))
        # print("Subgrid is ")
        # println(typeof(subgrid))
        # Count the number of elements in the subgrid
        cnt = 0
        for i in 1:length(grid)
            if issubindex(i, subgrid)
                cnt += 1
            end
        end
        @test cnt == length(subgrid)
    end

    g = subgrid(ScatteredGrid(rand(10)), Interval(0,.5))
    @test covering(g) ≈ Interval(0,.5)
    for x in g
        @test x ∈ covering(g)
    end

    for x in boundary(EquispacedGrid(100,-1,1)^2,UnitDisk())
        @test norm(x)≈1
    end

    for x in  boundary(subgrid(EquispacedGrid(100,-1,1)^2,UnitDisk()),UnitDisk())
        @test norm(x)≈1
    end
end
