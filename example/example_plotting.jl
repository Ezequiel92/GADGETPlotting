############################################################################################
# PLOTTING FUNCTIONS.
############################################################################################

scatterGridPlot(pos)
savefig(joinpath(BASE_OUT_PATH, "test_scatterGridPlot.png"))

densityMapPlot(pos, gas_mass, density, hsml, plane = "XY")
savefig(joinpath(BASE_OUT_PATH, "test_densityMapPlot_XY.png"))

densityMapPlot(pos, gas_mass, density, hsml, axes = true)
savefig(joinpath(BASE_OUT_PATH, "test_densityMapPlot_All.png"))

color = [:batlow, :bone, :CMRmap, :grayC, :inferno, :seaborn_rocket_gradient, :YlOrRd_9]

for c in color
    densityMapPlot(pos, gas_mass, density, hsml, color = c, axes = true)
    savefig(joinpath(BASE_OUT_PATH, "test_densityMapPlot_" * string(c) * ".png"))
end

starMapPlot(pos, plane = "XY", axes = true)
savefig(joinpath(BASE_OUT_PATH, "test_starMapPlot_XY.png"))

starMapPlot(pos, axes = true)
savefig(joinpath(BASE_OUT_PATH, "test_starMapPlot_All.png"))

gasStarEvolutionPlot(time_series, pos, SNAP_N)
savefig(joinpath(BASE_OUT_PATH, "test_gasStarEvolutionPlot.png"))

pgfplotsx()

CMDFPlot(star_mass, star_z, 1UnitfulAstro.Myr)
savefig(joinpath(BASE_OUT_PATH, "test_CMDFPlot.png"))

CMDFPlot([star_mass, star_mass], [star_z, star_z], 1UnitfulAstro.Myr, ["sim1" "sim2"])
savefig(joinpath(BASE_OUT_PATH, "test_compare_CMDFPlot.png"))

birthHistogramPlot(birth_pos, bins = 50)
savefig(joinpath(BASE_OUT_PATH, "test_birthHistogramPlot.png"))

timeSeriesPlot(time_series, mass_factor = 0, number_factor = 4)
savefig(joinpath(BASE_OUT_PATH, "test_timeSeriesPlot.png"))

scaleFactorSeriesPlot(time_series, mass_factor = 0, number_factor = 4)
savefig(joinpath(BASE_OUT_PATH, "test_scaleFactorSeriesPlot.png"))

redshiftSeriesPlot(time_series, mass_factor = 0, number_factor = 4)
savefig(joinpath(BASE_OUT_PATH, "test_redshiftSeriesPlot.png"))

compareSimulationsPlot(
    [time_series, time_series],
    "star_mass",
    "sfr",
    ["sim1" "sim2"];
    title = "SFMS relation",
    x_factor = 10,
    scale = (:identity, :log10),
)
savefig(joinpath(BASE_OUT_PATH, "test_compareSimulationsPlot.png"))

densityHistogramPlot(density, 1UnitfulAstro.Myr, factor = 10)
savefig(joinpath(BASE_OUT_PATH, "test_densityHistogramPlot.png"))

densityProfilePlot(pos, star_mass, 1UnitfulAstro.Myr, scale = :log10, bins = 50)
savefig(joinpath(BASE_OUT_PATH, "test_star_densityProfilePlot.png"))

densityProfilePlot(pos, gas_mass, 1UnitfulAstro.Myr, scale = :log10, bins = 50, factor = 6)
savefig(joinpath(BASE_OUT_PATH, "test_gas_densityProfilePlot.png"))

densityProfilePlot(
    [pos, pos],
    [gas_mass, gas_mass],
    1UnitfulAstro.Myr,
    ["sim_1" "sim_2"],
    bins = 50,
    scale = :log10,
    factor = 6,
)
savefig(joinpath(BASE_OUT_PATH, "test_compare_gas_densityProfilePlot.png"))

densityProfilePlot(
    [pos, pos],
    [star_mass, star_mass],
    1UnitfulAstro.Myr,
    ["sim_1" "sim_2"],
    bins = 50,
    scale = :log10,
)
savefig(joinpath(BASE_OUT_PATH, "test_compare_star_densityProfilePlot.png"))

metallicityProfilePlot(pos, gas_mass, gas_z, 1UnitfulAstro.Myr, scale = :log10, bins = 50)
savefig(joinpath(BASE_OUT_PATH, "test_gas_metallicityProfilePlot.png"))

metallicityProfilePlot(pos, star_mass, star_z, 1UnitfulAstro.Myr, scale = :log10, bins = 50)
savefig(joinpath(BASE_OUT_PATH, "test_star_metallicityProfilePlot.png"))

metallicityProfilePlot(
    [pos, pos],
    [gas_mass, gas_mass],
    [gas_z, gas_z],
    1UnitfulAstro.Myr,
    ["sim_1" "sim_2"],
    scale = :log10,
    bins = 50,
)
savefig(joinpath(BASE_OUT_PATH, "test_compare_gas_metallicityProfilePlot.png"))

metallicityProfilePlot(
    [pos, pos],
    [star_mass, star_mass],
    [star_z, star_z],
    1UnitfulAstro.Myr,
    ["sim_1" "sim_2"],
    scale = :log10,
    bins = 50,
)
savefig(joinpath(BASE_OUT_PATH, "test_compare_star_metallicityProfilePlot.png"))

massProfilePlot(pos, star_mass, 1UnitfulAstro.Myr, scale = :log10, bins = 50)
savefig(joinpath(BASE_OUT_PATH, "test_star_massProfilePlot.png"))

massProfilePlot(pos, gas_mass, 1UnitfulAstro.Myr, scale = :log10, bins = 50, factor = 10)
savefig(joinpath(BASE_OUT_PATH, "test_gas_massProfilePlot.png"))

massProfilePlot(
    [pos, pos],
    [gas_mass, gas_mass],
    1UnitfulAstro.Myr,
    ["sim_1" "sim_2"],
    scale = :log10,
    bins = 50,
    factor = 10,
)
savefig(joinpath(BASE_OUT_PATH, "test_compare_gas_massProfilePlot.png"))

massProfilePlot(
    [pos, pos],
    [star_mass, star_mass],
    1UnitfulAstro.Myr,
    ["sim_1" "sim_2"],
    scale = :log10,
    bins = 50,
)
savefig(joinpath(BASE_OUT_PATH, "test_compare_star_massProfilePlot.png"))

sfrTxtPlot(
    sfrtxt_data,
    1,
    [4, 6],
    title = "run_A_01",
    bins = 50,
    scale = (:identity, :log10),
)
savefig(joinpath(BASE_OUT_PATH, "test_sfrTxtPlot.png"))

temperatureHistogramPlot(temp_data, 1UnitfulAstro.Myr, bins = 30)
savefig(joinpath(BASE_OUT_PATH, "test_temperatureHistogramPlot.png"))

rhoTempPlot(temp_data, density, 1UnitfulAstro.Myr)
savefig(joinpath(BASE_OUT_PATH, "test_rhoTempPlot.png"))

KennicuttSchmidtPlot(
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
savefig(joinpath(BASE_OUT_PATH, "test_KennicuttSchmidtPlot.png"))