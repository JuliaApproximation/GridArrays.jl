
@testset "Plots" begin
    plot(FourierGrid(4))
    plot(FourierGrid(4)× FourierGrid(4) )
    plot(FourierGrid(4)× FourierGrid(4) × FourierGrid(4) )
    plot(FourierGrid(4)× FourierGrid(4) ,rand(4,4))
    plot(FourierGrid(4) ,rand(4))
end
