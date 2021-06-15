############################################################################################
# Test auxiliary functions
############################################################################################

##############
# Data
##############

pos_gas = pos["gas"]
pos_stars = pos["stars"]
m_gas = gas_mass["mass"]
dist_gas_2D_3D = sqrt.(pos_gas[1, :] .^ 2 .+ pos_gas[2, :] .^ 2 .+ pos_gas[3, :] .^ 2)
dist_gas_2D = sqrt.(pos_gas[1, :] .^ 2 + pos_gas[2, :] .^ 2)
dist_stars = sqrt.(pos_stars[1, :] .^ 2 + pos_stars[2, :] .^ 2)
box_size = ustrip(Float64, pos["unit"], BOX_SIZE)
max_z = findmax(star_z["Z"])
max_Z = max_z[1] / star_mass["mass"][max_z[2]] 
# To test comparison() and deep_comparison()
arr1 = [1.5, 9.6, 6.4]
arr2 = [1.5, 9.6, 6.9]
dict1 = Dict("a" => arr1, "b" => arr1)
dict2 = Dict("a" => arr1, "b" => arr2)

##############
# Functions
##############

relative_2D = GADGETPlotting.relative(plot(rand(100)), 0.5, 0.5)
relative_3D = GADGETPlotting.relative(surface(rand(100, 100)), 0.5, 0.5, 0.5)

# don't know how to test make_video() yet

smooth_w = GADGETPlotting.smooth_window([1:1000...], rand(1000), 50)

density_p = GADGETPlotting.density_profile(m_gas, dist_gas_2D_3D, box_size, 80)

metallicity_p = GADGETPlotting.metallicity_profile(
    m_gas, 
    dist_gas_2D_3D, 
    gas_z["Z"], 
    box_size, 
    80,
)

mass_p = GADGETPlotting.mass_profile(m_gas, dist_gas_2D_3D, box_size, 80)

cmdf = GADGETPlotting.compute_cmdf(star_mass["mass"], star_z["Z"], max_Z, 50)

ksl = GADGETPlotting.kennicutt_schmidt_law(
    m_gas,
    dist_gas_2D,
    temp_data["temperature"],
    star_mass["mass"],
    dist_stars,
    age_data["ages"],
    ustrip(Float64, temp_data["unit"], 3e4Unitful.K),
    ustrip(Float64, age_data["unit"], 200UnitfulAstro.Myr),
    ustrip(Float64, pos["unit"], BOX_SIZE),
    bins = 80,
)

# `FMOL` is not a block present in the examples of cosmological snapshots
if SIM_COSMO == 0
    quantities2D = GADGETPlotting.quantities_2D(
        m_gas,
        dist_gas_2D,
        temp_data["temperature"],
        star_mass["mass"],
        dist_stars,
        age_data["ages"],
        gas_mz["Z"],
        fmol,
        ustrip(Float64, temp_data["unit"], 3e4Unitful.K),
        ustrip(Float64, age_data["unit"], 200UnitfulAstro.Myr),	
        ustrip(Float64, pos["unit"], BOX_SIZE),
        bins = 80,
    )
end

format_err_1 = GADGETPlotting.format_error(69.42069, 0.038796)
format_err_2 = GADGETPlotting.format_error(69.42069, 0.018796)
format_err_3 = GADGETPlotting.format_error(69.42069, 0.0)
format_err_4 = GADGETPlotting.format_error(69.42069, 73.4)

pass = GADGETPlotting.pass_all(FIRST_SNAP, "gas")

energy_i = GADGETPlotting.energy_integrand(read_header(FIRST_SNAP), 1.0)

num_int_1 = GADGETPlotting.num_integrate(sin, 0, 3π)
num_int_2 = GADGETPlotting.num_integrate(x -> x^3 + 6 * x^2 + 9 * x + 2, 0, 4.69)
num_int_3 = GADGETPlotting.num_integrate(x -> exp(x^x), 0, 1)
num_int_4 = GADGETPlotting.num_integrate(x -> sqrt(sqrt(1 / (x + 1))), 0, 1)

rc = GADGETPlotting.center_of_mass(pos_gas, m_gas)

max_gas_dist = GADGETPlotting.max_length(pos_gas)

##############
# Testing 
##############

@testset "Auxiliary functions" begin

    temp_img = joinpath(@__DIR__, "test_img.png")
    gr()
    
    vline_plot = @test_nowarn GADGETPlotting.set_vertical_flags(
        ([4.0, 6.0], ["test_1", "test_2"]), 
        plot(1:10),
    )
    savefig(vline_plot, temp_img)

    jldopen(joinpath(BASE_DATA_PATH, "data_auxiliary.jld2"), "r") do file
        @test comparison(file["relative_2D"], relative_2D, rtol = 0.1)
        @test comparison(file["relative_3D"], relative_3D, rtol = 0.01)
        @test comparison(file["smooth_w"][1], smooth_w[1])
        @test deep_comparison(file["density_p"], density_p)
        @test deep_comparison(file["metallicity_p"], metallicity_p)
        @test deep_comparison(file["mass_p"], mass_p)
        @test deep_comparison(file["cmdf"], cmdf)
        if ksl === nothing
            @test file["ksl"] === ksl
        else
            @test comparison(file["ksl"]["RHO"], ksl["RHO"])
            @test comparison(file["ksl"]["SFR"], ksl["SFR"])
        end
        if SIM_COSMO == 0
            display(file["quantities2D"])
            display(quantities2D)
            @test deep_comparison(file["quantities2D"], quantities2D)
        end
        @test file["format_err_1"] == format_err_1
        @test file["format_err_2"] == format_err_2
        @test file["format_err_3"] == format_err_3
        @test file["format_err_4"] == format_err_4
        @test file["pass"] == pass
        @test file["energy_i"] ≈ energy_i
        @test file["num_int_1"] ≈ num_int_1
        @test file["num_int_2"] ≈ num_int_2
        @test file["num_int_3"] ≈ num_int_3
        @test file["num_int_4"] ≈ num_int_4
        @test comparison(file["rc"], rc)
        @test file["max_gas_dist"] ≈ max_gas_dist
        @test !comparison("test1", "test2")
        @test comparison(arr1, arr1)
        @test !deep_comparison(arr1, arr2)
        @test deep_comparison(arr1, arr1)
        @test !deep_comparison(dict1, dict2)
        @test deep_comparison(dict1, dict1)
        @test_reference joinpath(BASE_DATA_PATH, "vline_plot.png") load(temp_img)
    end

    @test_nowarn rm(temp_img)
	
end