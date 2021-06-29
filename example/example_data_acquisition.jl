############################################################################################
# Example data acquisition functions
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

if SIM_COSMO == 0

    jldsave(
        joinpath(BASE_OUT_PATH, "data_acquisition.jld2"); 
        snaps, time_series, pos, density, hsml, 
        gas_mass, dm_mass, star_mass, gas_z, star_z, 
        gas_mz, star_mz, temp_data, age_data, 
		birth_pos, sfr_txt_data, cpu_txt_data, 
        fmol, fatom,
    )

else

# Because of the GitHub 100MB per file limit, `pos`, `gas_mz`, `star_mz` and `cpu_txt_data`
# for are not included for cosmological simulations, and the rest of parameters are divided 
# in two files 

    jldsave(
        joinpath(BASE_OUT_PATH, "data_acquisition_01.jld2"); 
        snaps, time_series, density, hsml, gas_mass, dm_mass,
    )

    jldsave(
        joinpath(BASE_OUT_PATH, "data_acquisition_02.jld2"); 
        star_mass, gas_z, star_z, temp_data, 
        age_data, birth_pos, sfr_txt_data,
    )

end
