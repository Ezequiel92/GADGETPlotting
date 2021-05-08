push!(LOAD_PATH, "./src/")
using GADGETPlotting, GadgetIO, Test, ReferenceTests, Plots, JLD2, UnitfulAstro, Unitful

const BASE_SRC_PATH = joinpath(@__DIR__, "../example/example_data")
const FIRST_SNAP = joinpath(BASE_SRC_PATH, "snap_000")
const BASE_DATA_PATH = joinpath(@__DIR__, "../example/example_results")
const SNAP_NAME = "snap"
const BOX_SIZE = 200UnitfulAstro.kpc
const SIM_COSMO = 0
const SNAP_N = 21

@testset "GADGETPlotting tests" begin

    ########################################################################################
	# Test data acquisition functions
	########################################################################################

	snaps = GADGETPlotting.get_snapshot_path(SNAP_NAME, BASE_SRC_PATH)

	snap_files = snaps["snap_files"]
    sim_cosmo = SIM_COSMO
    snap_n = snap_files[SNAP_N]

	time_series = GADGETPlotting.get_time_evolution(snap_files; sim_cosmo)
	pos = GADGETPlotting.get_position(
		snap_n;
		sim_cosmo,
		box_size = BOX_SIZE,
		length_unit = UnitfulAstro.Mpc,
	)
	density = GADGETPlotting.get_density(snap_n; sim_cosmo)
	hsml = GADGETPlotting.get_hsml(snap_n; sim_cosmo)
	gas_mass = GADGETPlotting.get_mass(snap_n, "gas"; sim_cosmo)
	dm_mass = GADGETPlotting.get_mass(snap_n, "dark_matter"; sim_cosmo)
	star_mass = GADGETPlotting.get_mass(snap_n, "stars"; sim_cosmo)
	gas_z = GADGETPlotting.get_metallicity(snap_n, "gas"; sim_cosmo)
	star_z = GADGETPlotting.get_metallicity(snap_n, "stars"; sim_cosmo)
	temp_data = GADGETPlotting.get_temperature(snap_n; sim_cosmo)
	age_data = GADGETPlotting.get_age(
		snap_n,
		time_series["clock_time"][SNAP_N] * time_series["units"]["time"];
		sim_cosmo,
	)
	birth_pos = GADGETPlotting.get_birth_place(
		SNAP_N,
		snap_files,
		time_series["clock_time"],
		time_series["units"]["time"];
		sim_cosmo,
	)
	sfrtxt_data = GADGETPlotting.get_sfr_txt(BASE_SRC_PATH, FIRST_SNAP; sim_cosmo)
	
    jldopen(joinpath(BASE_DATA_PATH, "data_acquisition.jld2"), "r") do file
        @test file["snaps"]["numbers"] == snaps["numbers"]
        @test deep_comparison(file["time_series"], time_series)
		@test deep_comparison(file["pos"], pos)
        @test deep_comparison(file["density"], density)
        @test deep_comparison(file["hsml"], hsml)
        @test deep_comparison(file["gas_mass"], gas_mass)
        @test deep_comparison(file["dm_mass"], dm_mass)
        @test deep_comparison(file["star_mass"], star_mass)
        @test deep_comparison(file["gas_z"], gas_z)
        @test deep_comparison(file["star_z"], star_z)
        @test deep_comparison(file["temp_data"], temp_data)
        @test deep_comparison(file["age_data"], age_data)
        @test deep_comparison(file["birth_pos"], birth_pos)
        @test deep_comparison(file["sfrtxt_data"], sfrtxt_data)
    end

    ########################################################################################
    # Test plotting functions
    ########################################################################################

    temp_img = joinpath(@__DIR__, "test_img.png")

    @test_nowarn fig = scatter_grid_plot(pos)
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "scatter_grid_plot.png") load(temp_img)

    @test_nowarn fig = density_map_plot(pos, gas_mass, density, hsml)
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "density_map_plot.png") load(temp_img)

    @test_nowarn fig = star_map_plot(pos)
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "star_map_plot_All.png") load(temp_img)

    @test_nowarn fig = gas_star_evolution_plot(SNAP_N, time_series, pos)
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "gas_star_evolution_plot.png") load(temp_img)

    @test_nowarn fig = cmdf_plot(star_mass, star_z, 1UnitfulAstro.Myr)
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "cmdf_plot.png") load(temp_img)

    @test_nowarn fig = cmdf_plot(
        [star_mass, star_mass], 
        [star_z, star_z], 
        1UnitfulAstro.Myr, 
        ["sim1" "sim2"],
    )
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "compare_cmdf_plot.png") load(temp_img)

    @test_nowarn fig = birth_histogram_plot(birth_pos, bins = 50)
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "birth_histogram_plot.png") load(temp_img)

    @test_nowarn fig = time_series_plot(time_series, mass_factor = 0, number_factor = 4)
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "time_series_plot.png") load(temp_img)

    @test_nowarn fig = scale_factor_series_plot(
        time_series, 
        mass_factor = 0, 
        number_factor = 4,
    )
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "scale_factor_series_plot.png") load(temp_img)

    @test_nowarn fig = redshift_series_plot(time_series, mass_factor = 0, number_factor = 4)
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "redshift_series_plot.png") load(temp_img)

    @test_nowarn fig = compare_simulations_plot(
        [time_series, time_series],
        "star_mass",
        "sfr",
        ["sim1" "sim2"];
        title = "SFMS relation",
        x_factor = 10,
        scale = [:identity, :log10],
    )
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "compare_simulations_plot.png") load(temp_img)

    @test_nowarn fig = density_histogram_plot(density, 1UnitfulAstro.Myr, factor = 10)
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "density_histogram_plot.png") load(temp_img)

    @test_nowarn fig = density_profile_plot(
        pos, 
        gas_mass, 
        1UnitfulAstro.Myr, 
        scale = :log10, 
        bins = 50, 
        factor = 6,
    )
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "density_profile_plot.png") load(temp_img)
    
    @test_nowarn fig = density_profile_plot(
        [pos, pos],
        [gas_mass, gas_mass],
        1UnitfulAstro.Myr,
        ["sim_1" "sim_2"],
        bins = 50,
        scale = :log10,
        factor = 6,
    )
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "compare_density_profile_plot.png") load(temp_img)

    @test_nowarn fig = metallicity_profile_plot(
        pos, 
        gas_mass, 
        gas_z, 
        1UnitfulAstro.Myr, 
        scale = :log10, 
        bins = 50,
    )
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "metallicity_profile_plot.png") load(temp_img)

    @test_nowarn fig = metallicity_profile_plot(
        [pos, pos],
        [gas_mass, gas_mass],
        [gas_z, gas_z],
        1UnitfulAstro.Myr,
        ["sim_1" "sim_2"],
        scale = :log10,
        bins = 50,
    )
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "compare_metallicity_profile_plot.png") load(temp_img)

    @test_nowarn fig = mass_profile_plot(
        pos, 
        gas_mass, 
        1UnitfulAstro.Myr, 
        scale = :log10, 
        bins = 50, 
        factor = 10,
    )
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "mass_profile_plot.png") load(temp_img)

    @test_nowarn fig = mass_profile_plot(
        [pos, pos],
        [gas_mass, gas_mass],
        1UnitfulAstro.Myr,
        ["sim_1" "sim_2"],
        scale = :log10,
        bins = 50,
        factor = 10,
    )
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "compare_mass_profile_plot.png") load(temp_img)

    @test_nowarn fig = sfr_txt_plot(
        sfrtxt_data,
        1,
        [4, 6],
        title = "run_A_01",
        bins = 50,
        scale = (:identity, :log10),
    )
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "compare_columns_sfr_txt_plot.png") load(temp_img)

    @test_nowarn fig = sfr_txt_plot(
        [sfrtxt_data, sfrtxt_data],
        1,
        6,
        ["sim_1" "sim_2"],
        title = "Column 6 vs 1",
        bins = 50,
        scale = (:identity, :log10),
    )
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "compare_sims_sfr_txt_plot.png") load(temp_img)

    @test_nowarn fig = temperature_histogram_plot(temp_data, 1UnitfulAstro.Myr, bins = 30)
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "temperature_histogram_plot.png") load(temp_img)

    @test_nowarn fig = rho_temp_plot(temp_data, density, 1UnitfulAstro.Myr)
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "rho_temp_plot.png") load(temp_img)

    @test_nowarn fig = kennicutt_schmidt_plot(
        gas_mass,
        temp_data,
        star_mass,
        age_data,
        pos,
        3.0e4Unitful.K,
        20UnitfulAstro.Myr,
        BOX_SIZE,
        1UnitfulAstro.Myr,
        bins = 80,
        error_formating = "conf_interval",
    )
    # Base.invokelatest(savefig, fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "kennicutt_schmidt_plot.png") load(temp_img)
    
    # @test_nowarn rm(temp_img)

    ########################################################################################
	# Test auxiliary functions
    ########################################################################################

    relative_2D = GADGETPlotting.relative(plot(rand(100)), 0.5, 0.5)
    relative_3D = GADGETPlotting.relative(surface(rand(100, 100)), 0.5, 0.5, 0.5)

    smooth_w = GADGETPlotting.smooth_window([1:1000...], rand(1000), 50)

    positions = pos["gas"]
    distances = sqrt.(positions[1, :] .^ 2 .+ positions[2, :] .^ 2 .+ positions[3, :] .^ 2)
    box_size = ustrip(Float64, pos["unit"], BOX_SIZE)
    density_p = GADGETPlotting.density_profile(gas_mass["mass"], distances, box_size, 80)

    metallicity_p = GADGETPlotting.metallicity_profile(
        gas_mass["mass"], 
        distances, 
        gas_z["Z"], 
        box_size, 
        80,
    )

    mass_p = GADGETPlotting.mass_profile(gas_mass["mass"], distances, box_size, 80)

    max_z = findmax(star_z["Z"])
    max_Z = max_z[1] / star_mass["mass"][max_z[2]]
    cmdf = GADGETPlotting.compute_cmdf(star_mass["mass"], star_z["Z"], max_Z, 50)

    pos_gas = pos["gas"]
    dist_gas = sqrt.(pos_gas[1, :] .^ 2 + pos_gas[2, :] .^ 2)
    pos_stars = pos["stars"]
    dist_stars = sqrt.(pos_stars[1, :] .^ 2 + pos_stars[2, :] .^ 2)
    ksl = GADGETPlotting.kennicutt_schmidt_law(
        gas_mass["mass"],
        dist_gas,
        temp_data["temperature"],
        star_mass["mass"],
        dist_stars,
        age_data["ages"],
        ustrip(Float64, temp_data["unit"], 3e4Unitful.K),
        ustrip(Float64, age_data["unit"], 200UnitfulAstro.Myr),
        ustrip(Float64, pos["unit"], BOX_SIZE),
        bins = 80,
    )

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

    # To test comparison() and deep_comparison()
    arr1 = [1.5, 9.6, 6.4]
    arr2 = [1.5, 9.6, 6.9]
    dict1 = Dict("a" => arr1, "b" => arr1)
    dict2 = Dict("a" => arr1, "b" => arr2)

    jldopen(joinpath(BASE_DATA_PATH, "data_auxiliary.jld2"), "r") do file
        @test comparison(file["relative_2D"], relative_2D, rtol = 0.1)
        @test comparison(file["relative_3D"], relative_3D, rtol = 0.1)
        @test comparison(file["smooth_w"][1], smooth_w[1])
        @test deep_comparison(file["density_p"], density_p)
        @test deep_comparison(file["metallicity_p"], metallicity_p)
        @test deep_comparison(file["mass_p"], mass_p)
        @test deep_comparison(file["cmdf"], cmdf)
        @test comparison(file["ksl"]["RHO"], ksl["RHO"])
        @test comparison(file["ksl"]["SFR"], ksl["SFR"])
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
        @test !comparison("test1", "test2")
        @test comparison(arr1, arr1)
        @test !deep_comparison(arr1, arr2)
        @test deep_comparison(arr1, arr1)
        @test !deep_comparison(dict1, dict2)
        @test deep_comparison(dict1, dict1)
    end
	
end