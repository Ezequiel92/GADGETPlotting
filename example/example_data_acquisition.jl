############################################################################################
# Example data acquisition functions
############################################################################################

##############
## Functions
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

# The cpu.txt of cosmological simulations is too heavy for GitHub, and `get_cpu_txt`
# doesn't distinguish among different types of simulations.
if SIM_COSMO == 0
    cpu_txt_data = GADGETPlotting.get_cpu_txt(BASE_SRC_PATH, ["i/o", "hotngbs", "density"])
end

##############
## Testing 
##############

# No `pos` and two files for cosmological simulations, 
# because of the GitHub 100MB limit per file

if SIM_COSMO == 0

    jldsave(
        joinpath(BASE_OUT_PATH, "data_acquisition.jld2"); 
        snaps, time_series, pos, density, hsml, 
        gas_mass, dm_mass, star_mass, gas_z, star_z, 
        temp_data, age_data, birth_pos, sfr_txt_data,
        cpu_txt_data,
    )

else

    jldsave(
        joinpath(BASE_OUT_PATH, "data_acquisition_01.jld2"); 
        snaps, time_series, density, hsml, gas_mass, 
    )

    jldsave(
        joinpath(BASE_OUT_PATH, "data_acquisition_02.jld2"); 
        dm_mass, star_mass, gas_z, star_z, temp_data, 
        age_data, birth_pos, sfr_txt_data,
    )

end
