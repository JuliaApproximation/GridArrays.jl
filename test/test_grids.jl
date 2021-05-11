
interval_grids = [EquispacedGrid,
    PeriodicEquispacedGrid, MidpointEquispacedGrid,
    UnitEquispacedGrid, UnitMidpointEquispacedGrid,
    UnitPeriodicEquispacedGrid,
    ChebyshevNodes, ChebyshevExtremae,
    ChebyshevUNodes, LegendreNodes]

function test_grids(T)
    ## Equispaced grids
    len = 21
    a = -T(1.2)
    b = T(3.5)
    g1 = EquispacedGrid(len, Interval(a, b))
    g2 = PeriodicEquispacedGrid(len, Interval(a-1, b+1))
    g3 = g1 × g2
    g4 = g1 × g3

    test_generic_grid(g1)
    test_generic_grid(g2)
    test_generic_grid(g3)
    test_generic_grid(g4)

    # Test a subgrid
    g5 = g1[10:20]
    @test g5[1] == g1[10]
    @test g5[11] == g1[20]
    @test length(g5) == 20-10+1
    test_generic_grid(g5)
    g6 = g1[10:2:20]
    @test g6[2] == g1[12]
    @test length(g6) == 6

    g = EquispacedGrid(len, a, b)
    idx = 5
    @test g[idx] ≈ a + (idx-1) * (b-a)/(len-1)
    @test g[len] ≈ b
    @test_throws BoundsError g[len+1] == b

    # Test iterations
    (l,s) = grid_iterator1(g)
    @assert s ≈ len * (a+b)/2

    ## Periodic equispaced grids
    len = 120
    a = -T(1.2)
    b = T(3.5)
    g = PeriodicEquispacedGrid(len, a, b)

    idx = 5
    @test g[idx] ≈ a + (idx-1) * (b-a)/len
    @test g[len] ≈ b - step(g)
    @test_throws BoundsError g[len+1] == b

    (l,s) = grid_iterator1(g)
    @test s ≈ (len-1)*(a+b)/2 + a

    (l,s) = grid_iterator2(g)
    @test s ≈ (len-1)*(a+b)/2 + a

    ## Tensor product grids
    len = 11
    g1 = PeriodicEquispacedGrid(len, -one(T), one(T))
    g = g1^2
    @test isperiodic(g1)
    @test isperiodic(g)
    @test size(g) == (length(g1),length(g1))
    @test productgrid(g1) ==g1

    g2 = EquispacedGrid(len, -one(T), one(T))
    @test !isperiodic(g2)
    g = g1 × g2
    @test !isperiodic(g)
    @test length(g) == length(g1) * length(g2)
    @test size(g) == (length(g1),length(g2))
    @test size(g,1) == length(g1)

    @test component(g, 1) == g1
    @test component(g, 2) == g2
    @test covering(g) ≈ covering(g1)×covering(g2)



    idx1 = 5
    idx2 = 9
    x1 = g1[idx1]
    x2 = g2[idx2]
    x = g[idx1,idx2]
    @test x[1] ≈ x1
    @test x[2] ≈ x2

    (l,s) = grid_iterator1(g)
    @test s ≈ -len

    (l,s) = grid_iterator2(g)
    @test s ≈ -len

    # Test a tensor of a tensor
    g3 = g × g2
    idx1 = 5
    idx2 = 7
    idx3 = 4
    x = g3[idx1,idx2,idx3]
    x1 = g1[idx1]
    x2 = g2[idx2]
    x3 = g2[idx3]
    @test x[1] ≈ x1
    @test x[2] ≈ x2
    @test x[3] ≈ x3

    # Test a mapped grid
    m = mapto(T(0)..T(1), T(2)..T(3))
    # Make a MappedGrid by hand because map_grid would simplify
    mg1 = MappedGrid(m, PeriodicEquispacedGrid(30, T(0), T(1)))
    test_generic_grid(mg1)
    # Does map_grid simplify?
    mg2 = map_grid(m, PeriodicEquispacedGrid(30, T(0), T(1)))
    @test typeof(mg2) <: PeriodicEquispacedGrid
    @test infimum(covering(mg2)) ≈ T(2)
    @test supremum(covering(mg2)) ≈ T(3)

    # Apply a second map and check whether everything simplified
    m2 = mapto(T(2)..T(3), T(4)..T(5))
    mg3 = map_grid(m2, mg1)
    @test infimum(covering(mg3)) ≈ T(4)
    @test supremum(covering(mg3)) ≈ T(5)
    @test mg3 isa PeriodicEquispacedGrid

    # Scattered grid
    pts = rand(T, 10)
    sg = ScatteredGrid(pts)
    test_generic_grid(sg)
end

function test_laguerre(T)
    grid = LaguerreNodes(10,rand(T))
    @test infimum(covering(LaguerreNodes(0.,rand(10)))) == 0
    @test supremum(covering(LaguerreNodes(0.,rand(10)))) == Inf
    test_generic_grid(grid)
end

function test_hermite(T)
    grid = HermiteNodes(10)
    @test DomainSets.FullSpace{T}()== covering(HermiteNodes(rand(T,10)))
    test_generic_grid(grid)
end

function test_jacobi(T)
    grid = JacobiNodes(10,rand(T),rand(T))
    @test covering(JacobiNodes(T(0),T(0),rand(T,10))) == ChebyshevInterval{T}()
    test_generic_grid(grid)
    @test JacobiNodes(10,zero(T),zero(T)) ≈ LegendreNodes(10)
    @test JacobiNodes(10,T(1//2),T(1//2)) ≈ ChebyshevUNodes(10)
    @test JacobiNodes(10,T(-1//2),T(-1//2)) ≈ ChebyshevTNodes(10)
end



for T in types
    delimit(string(T))
    for GRID in interval_grids
        @testset "$(string(GRID))" begin
            if GRID <: GridArrays.AbstractEquispacedRangeGrid
                g = GRID(10,UnitInterval{T}())
            else
                g = GRID(10)
            end
            test_interval_grid(g)
        end
    end

    for grid in (JacobiNodes(10,rand(T),rand(T)),)
        @testset "JacobiNodes" begin
            test_interval_grid(grid)
        end
    end

    @testset "HermiteNodes" begin
        test_hermite(T)
    end

    @testset "Laguerre" begin
        test_laguerre(T)
    end

    @testset "Jacobi" begin
        test_jacobi(T)
    end

    @testset "Specific grid tests" begin
        test_grids(T)
    end
end
