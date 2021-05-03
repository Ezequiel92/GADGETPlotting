push!(LOAD_PATH, "./src/")
using GADGETPlotting, GadgetIO, Test, ReferenceTests, Plots

@testset "Unit Conversion" begin
	
	# Test auxiliary functions
	
	fig2D = plot(rand(100))
	@test all(.≈(
        GADGETPlotting.relative(fig2D, 0.5, 0.5), 
        (50.5, 0.5), 
        atol = 0.1,
    ))
	fig3D = surface(rand(100, 100))
	@test all(.≈(
        GADGETPlotting.relative(fig3D, 0.5, 0.5, 0.5), 
        (50.5, 50.5, 0.5), 
        atol = 0.1,
    ))
	
	x_data = [1:1000...]
	y_data = rand(1000)
    x_smooth = GADGETPlotting.smoothWindow(x_data, y_data, 50)[1]
	@test x_smooth == [10.5, 30.5, 50.5, 70.5, 90.5, 110.5, 130.5, 150.5, 170.5, 190.5, 210.5, 230.5, 250.5, 270.5, 290.5, 310.5, 330.5, 350.5, 370.5, 390.5, 410.5, 430.5, 450.5, 470.5, 490.5, 510.5, 530.5, 550.5, 570.5, 590.5, 610.5, 630.5, 650.5, 670.5, 690.5, 710.5, 730.5, 750.5, 770.5, 790.5, 810.5, 830.5, 850.5, 870.5, 890.5, 910.5, 930.5, 950.5, 970.5, 990.0]
	
	@test GADGETPlotting.format_error(69.42069, 0.038796) == (69.42, 0.04)
	@test GADGETPlotting.format_error(69.42069, 0.018796) == (69.421, 0.019)
	@test GADGETPlotting.format_error(69.42069, 0.0) == (69.42069, 0.0)
	@test GADGETPlotting.format_error(69.42069, 73.4) == (0.0, 70.0)
	
    nums = GADGETPlotting.pass_all(
        joinpath(@__DIR__, "../example/test_data/snap_000"), 
        "gas",
    )
	@test length(nums) == 17256
    @test sum(nums) == 148893396
	
	header = read_header(joinpath(@__DIR__, "../example/test_data/snap_000"))
	@test GADGETPlotting.energy_integrand(header, 1.0) ≈ 9.7846401
	
	@test GADGETPlotting.num_integrate(sin, 0, 3π) ≈ 1.999629876136
	@test GADGETPlotting.num_integrate(x -> x^3 + 6 * x^2 + 9 * x + 2, 0, 4.69) ≈ 438.9004836
	@test GADGETPlotting.num_integrate(x -> exp(x^x), 0, 1) ≈ 2.1975912
	
end