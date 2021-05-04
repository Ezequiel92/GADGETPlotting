push!(LOAD_PATH, "./src/")
using GADGETPlotting, GadgetIO, Test, ReferenceTests, Plots, JLD2, UnitfulAstro, Unitful

const BASE_SRC_PATH = joinpath(@__DIR__, "../example/example_data")
const FIRST_SNAP = joinpath(BASE_SRC_PATH, "snap_000")
const BASE_DATA_PATH = joinpath(@__DIR__, "../example/example_results")
const SNAP_NAME = "snap"
const BOX_SIZE = 200UnitfulAstro.kpc
const SIM_COSMO = 0
const SNAP_N = 21

function vec_compare(x, y; rtol::Float64 = 1e-5)::Bool
    return all(.≈(x, y; rtol))
end

function deep_comparison(x::Dict, y::Dict)::Bool

    equal = Bool[]
    for key in keys(x)
        if typeof(x[key]) <: AbstractArray{<:Real}
            push!(equal, vec_compare(x[key], y[key]))
        else
            push!(equal, x[key] == y[key])
        end
    end

    return all(equal)
end

@testset "Unit Conversion" begin

    ########################################################################################
	# Test data acquisition
	########################################################################################

	snaps = GADGETPlotting.getSnapshotPaths(SNAP_NAME, BASE_SRC_PATH)

	snap_files = snaps["snap_files"]
    sim_cosmo = SIM_COSMO
    snap_n = snap_files[SNAP_N]

	time_series = GADGETPlotting.timeSeriesData(snap_files; sim_cosmo)

	pos = GADGETPlotting.positionData(
		snap_n;
		sim_cosmo,
		box_size = BOX_SIZE,
		length_unit = UnitfulAstro.Mpc,
	)

	density = GADGETPlotting.densityData(snap_n; sim_cosmo)
	hsml = GADGETPlotting.hsmlData(snap_n; sim_cosmo)
	gas_mass = GADGETPlotting.massData(snap_n, "gas"; sim_cosmo)
	dm_mass = GADGETPlotting.massData(snap_n, "dark_matter"; sim_cosmo)
	star_mass = GADGETPlotting.massData(snap_n, "stars"; sim_cosmo)
	gas_z = GADGETPlotting.zData(snap_n, "gas"; sim_cosmo)
	star_z = GADGETPlotting.zData(snap_n, "stars"; sim_cosmo)
	temp_data = GADGETPlotting.tempData(snap_n; sim_cosmo)
	age_data = GADGETPlotting.ageData(
		snap_n,
		time_series["clock_time"][SNAP_N] * time_series["units"]["time"];
		sim_cosmo,
	)
	birth_pos = GADGETPlotting.birthPlace(
		SNAP_N,
		snap_files,
		time_series["clock_time"],
		time_series["units"]["time"];
		sim_cosmo,
	)
	sfrtxt_data = GADGETPlotting.sfrTxtData(BASE_SRC_PATH, FIRST_SNAP; sim_cosmo)
	
    jldopen(joinpath(BASE_DATA_PATH, "data_acquisition.jld2"), "r") do file
        @test file["snaps"]["numbers"] == snaps["numbers"]
        @test deep_comparison(file["time_series"], time_series)
        @test vec_compare(file["pos"]["gas"], pos["gas"])
        @test vec_compare(file["pos"]["dark_matter"], pos["dark_matter"])
        @test vec_compare(file["pos"]["stars"], pos["stars"])
        @test file["pos"]["box_size"] == pos["box_size"]
        @test file["pos"]["periodic"] == pos["periodic"]
        @test file["pos"]["unit"] == pos["unit"]
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
  
	# Test auxiliary functions

    relative_2D = GADGETPlotting.relative(plot(rand(100)), 0.5, 0.5)
    relative_3D = GADGETPlotting.relative(surface(rand(100, 100)), 0.5, 0.5, 0.5)

    smooth_w = GADGETPlotting.smoothWindow([1:1000...], rand(1000), 50)

    positions = pos["gas"]
    distances = sqrt.(positions[1, :] .^ 2 .+ positions[2, :] .^ 2 .+ positions[3, :] .^ 2)
    box_size = ustrip(Float64, pos["unit"], BOX_SIZE)
    density_p = GADGETPlotting.densityProfile(gas_mass["mass"], distances, box_size, 80)

    metallicity_p = GADGETPlotting.metallicityProfile(
        gas_mass["mass"], 
        distances, 
        gas_z["Z"], 
        box_size, 
        80,
    )

    mass_p = GADGETPlotting.massProfile(gas_mass["mass"], distances, box_size, 80)

    max_z = findmax(star_z["Z"])
    max_Z = max_z[1] / star_mass["mass"][max_z[2]]
    cmdf = GADGETPlotting.CMDF(star_mass["mass"], star_z["Z"], max_Z, 50)

    pos_gas = pos["gas"]
    dist_gas = sqrt.(pos_gas[1, :] .^ 2 + pos_gas[2, :] .^ 2)
    pos_stars = pos["stars"]
    dist_stars = sqrt.(pos_stars[1, :] .^ 2 + pos_stars[2, :] .^ 2)
    ksl = GADGETPlotting.KennicuttSchmidtLaw(
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

    jldopen(joinpath(BASE_DATA_PATH, "data_auxiliary.jld2"), "r") do file
        @test vec_compare(file["relative_2D"], relative_2D, rtol = 0.1)
        @test vec_compare(file["relative_3D"], relative_3D, rtol = 0.1)
        @test vec_compare(file["smooth_w"][1], smooth_w[1])
        @test vec_compare(file["density_p"][1], density_p[1])
        @test vec_compare(file["density_p"][2], density_p[2])
        @test vec_compare(file["metallicity_p"][1], metallicity_p[1])
        @test vec_compare(file["metallicity_p"][2], metallicity_p[2])
        @test vec_compare(file["mass_p"][1], mass_p[1])
        @test vec_compare(file["mass_p"][2], mass_p[2])
        @test vec_compare(file["cmdf"][1], cmdf[1])
        @test vec_compare(file["cmdf"][2], cmdf[2])
        @test vec_compare(file["ksl"]["RHO"], ksl["RHO"])
        @test vec_compare(file["ksl"]["SFR"], ksl["SFR"])
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
    end
	
end