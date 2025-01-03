############################################################################################
# Test data acquisition functions
############################################################################################

##############
# Functions
##############

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

gas_mz = GADGETPlotting.get_metal_mass(snap_n, "gas"; sim_cosmo)

star_mz = GADGETPlotting.get_metal_mass(snap_n, "stars"; sim_cosmo)

temp_data = GADGETPlotting.get_temperature(snap_n; sim_cosmo)

age_data = GADGETPlotting.get_age(
	snap_n,
	time_series["clock_time"][SNAP_N] * time_series["units"]["time"];
	sim_cosmo,
    snap_0 = FIRST_SNAP,
)

birth_pos = GADGETPlotting.get_birth_place(
	SNAP_N,
	snap_files,
	time_series["clock_time"],
	time_series["units"]["time"];
	sim_cosmo,
)

sfr_txt_data = GADGETPlotting.get_sfr_txt(BASE_SRC_PATH, FIRST_SNAP; sim_cosmo)

if SIM_COSMO == 0
    
    # The file cpu.txt from a cosmological simulation is too big for GitHub, and `get_cpu_txt`
    # doesn't distinguish among different types of simulations, so there is no need to test it
    # when `SIM_COSMO` = 1
	cpu_txt_data = GADGETPlotting.get_cpu_txt(BASE_SRC_PATH, ["i/o", "hotngbs", "density"])

    # `FMOL` is not a block present in the examples of cosmological snapshots
    fmol = GADGETPlotting.get_fmol(snap_n; sim_cosmo)

    # `NH` is not a block present in the examples of cosmological snapshots
    fatom = GADGETPlotting.get_fatom(snap_n; sim_cosmo)

    quantities2D = GADGETPlotting.quantities_2D(
        gas_mass["mass"],
        sqrt.(pos["gas"][1, :] .^ 2 + pos["gas"][2, :] .^ 2),
        temp_data["temperature"],
        star_mass["mass"],
        sqrt.(pos["stars"][1, :] .^ 2 + pos["stars"][2, :] .^ 2),
        age_data["ages"],
        gas_mz["Z"],
        fmol,
        ustrip(Float64, temp_data["unit"], 3e4Unitful.K),
        ustrip(Float64, age_data["unit"], 200UnitfulAstro.Myr),	
        ustrip(Float64, pos["unit"], BOX_SIZE),
        bins = 80,
    )
end

##############
# Testing 
##############

@testset "Data acquisition functions" begin 

    if SIM_COSMO == 0
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
			@test deep_comparison(file["gas_mz"], gas_mz)
            @test deep_comparison(file["star_mz"], star_mz)
            @test deep_comparison(file["temp_data"], temp_data)
            @test deep_comparison(file["age_data"], age_data)
            @test deep_comparison(file["birth_pos"], birth_pos)
            @test deep_comparison(file["sfr_txt_data"], sfr_txt_data)
            @test deep_comparison(file["cpu_txt_data"], cpu_txt_data)
            @test comparison(file["fmol"], fmol)
            @test comparison(file["fatom"], fatom)
        end
    else

        # Because of the GitHub 100MB per file limit, `pos`, `gas_mz`, `star_mz` and `cpu_txt_data`
        # for are not included for cosmological simulations, and the rest of parameters are divided 
        # in two files 
        
        jldopen(joinpath(BASE_DATA_PATH, "data_acquisition_01.jld2"), "r") do file
            @test file["snaps"]["numbers"] == snaps["numbers"]
            @test deep_comparison(file["time_series"], time_series)
            @test deep_comparison(file["density"], density)
            @test deep_comparison(file["hsml"], hsml)
            @test deep_comparison(file["gas_mass"], gas_mass)
            @test deep_comparison(file["dm_mass"], dm_mass)
        end
        jldopen(joinpath(BASE_DATA_PATH, "data_acquisition_02.jld2"), "r") do file
            
            @test deep_comparison(file["star_mass"], star_mass)
            @test deep_comparison(file["gas_z"], gas_z)
            @test deep_comparison(file["star_z"], star_z)
            @test deep_comparison(file["temp_data"], temp_data)
            @test deep_comparison(file["age_data"], age_data)
            @test deep_comparison(file["birth_pos"], birth_pos)
            @test deep_comparison(file["sfr_txt_data"], sfr_txt_data)
        end
    end
	
end