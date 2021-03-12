############################################################################################
# DATA ACQUISITION FUNCTIONS.
############################################################################################

snaps = GADGETPlotting.getSnapshotPaths(SNAP_NAME, BASE_SRC_PATH)
snap_files = snaps["snap_files"]
display(snaps)
println()

time_series = GADGETPlotting.timeSeriesData(snap_files, sim_cosmo = SIM_COSMO)
display(time_series)
println()

pos = GADGETPlotting.positionData(
    snap_files[SNAP_N],
    sim_cosmo = SIM_COSMO,
    box_size = BOX_SIZE,
    length_unit = UnitfulAstro.Mpc,
)
display(pos)
println()

density = GADGETPlotting.densityData(snap_files[SNAP_N], sim_cosmo = SIM_COSMO)
display(density)
println()

hsml = GADGETPlotting.hsmlData(snap_files[SNAP_N], sim_cosmo = SIM_COSMO)
display(hsml)
println()

gas_mass = GADGETPlotting.massData(snap_files[SNAP_N], "gas", sim_cosmo = SIM_COSMO)
display(gas_mass)
println()

dm_mass = GADGETPlotting.massData(snap_files[SNAP_N], "dark_matter", sim_cosmo = SIM_COSMO)
display(dm_mass)
println()

star_mass = GADGETPlotting.massData(snap_files[SNAP_N], "stars", sim_cosmo = SIM_COSMO)
display(star_mass)
println()

gas_z = GADGETPlotting.zData(snap_files[SNAP_N], "gas", sim_cosmo = SIM_COSMO)
display(gas_z)
println()

star_z = GADGETPlotting.zData(snap_files[SNAP_N], "stars", sim_cosmo = SIM_COSMO)
display(star_z)
println()

temp_data = GADGETPlotting.tempData(snap_files[SNAP_N], sim_cosmo = SIM_COSMO)
display(temp_data)
println()

age_data = GADGETPlotting.ageData(
    snap_files[SNAP_N],
    time_series["clock_time"][SNAP_N] * time_series["units"]["time"],
    sim_cosmo = SIM_COSMO,
)
display(age_data)
println()

birth_pos = GADGETPlotting.birthPlace(
    SNAP_N,
    snap_files,
    time_series["clock_time"],
    time_series["units"]["time"],
    sim_cosmo = SIM_COSMO,
)
display(birth_pos)
println()

sfrtxt_data = GADGETPlotting.sfrTxtData(BASE_SRC_PATH, SNAP_NAME, sim_cosmo = SIM_COSMO)
display(sfrtxt_data)
println()
