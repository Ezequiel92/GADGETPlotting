############################################################################################
# DATA ACQUISITION FUNCTIONS
############################################################################################

snaps = GADGETPlotting.getSnapshotPaths(SNAP_NAME, BASE_SRC_PATH)
snap_files = snaps["snap_files"]

time_series = GADGETPlotting.timeSeriesData(snap_files, sim_cosmo = SIM_COSMO)

pos = GADGETPlotting.positionData(
    snap_files[SNAP_N],
    sim_cosmo = SIM_COSMO,
    box_size = BOX_SIZE,
    length_unit = UnitfulAstro.Mpc,
)

density = GADGETPlotting.densityData(snap_files[SNAP_N], sim_cosmo = SIM_COSMO)

hsml = GADGETPlotting.hsmlData(snap_files[SNAP_N], sim_cosmo = SIM_COSMO)

gas_mass = GADGETPlotting.massData(snap_files[SNAP_N], "gas", sim_cosmo = SIM_COSMO)

dm_mass = GADGETPlotting.massData(snap_files[SNAP_N], "dark_matter", sim_cosmo = SIM_COSMO)

star_mass = GADGETPlotting.massData(snap_files[SNAP_N], "stars", sim_cosmo = SIM_COSMO)

gas_z = GADGETPlotting.zData(snap_files[SNAP_N], "gas", sim_cosmo = SIM_COSMO)

star_z = GADGETPlotting.zData(snap_files[SNAP_N], "stars", sim_cosmo = SIM_COSMO)

temp_data = GADGETPlotting.tempData(snap_files[SNAP_N], sim_cosmo = SIM_COSMO)

age_data = GADGETPlotting.ageData(
    snap_files[SNAP_N],
    time_series["clock_time"][SNAP_N] * time_series["units"]["time"],
    sim_cosmo = SIM_COSMO,
)

birth_pos = GADGETPlotting.birthPlace(
    SNAP_N,
    snap_files,
    time_series["clock_time"],
    time_series["units"]["time"],
    sim_cosmo = SIM_COSMO,
)

sfrtxt_data = GADGETPlotting.sfrTxtData(
    BASE_SRC_PATH, 
    joinpath(BASE_SRC_PATH, SNAP_NAME * "_000"), 
    sim_cosmo = SIM_COSMO
)

# Save result for tests
jldsave(
    joinpath(BASE_OUT_PATH, "data_acquisition.jld2"); 
    snaps, time_series, pos, 
    density, hsml, gas_mass, 
    dm_mass, star_mass, gas_z, 
    star_z, temp_data, age_data, 
    birth_pos, sfrtxt_data,
)
