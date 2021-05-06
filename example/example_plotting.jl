############################################################################################
# PLOTTING FUNCTIONS
############################################################################################

figure = scatterGridPlot(pos)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "scatterGridPlot.svg"))

figure = densityMapPlot(pos, gas_mass, density, hsml)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "densityMapPlot.png"))

figure = starMapPlot(pos)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "starMapPlot.png"))

figure = gasStarEvolutionPlot(SNAP_N, time_series, pos)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "gasStarEvolutionPlot.png"))

figure = CMDFPlot(star_mass, star_z, 1UnitfulAstro.Myr)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "CMDFPlot.png"))

figure = CMDFPlot(
    [star_mass, star_mass], 
    [star_z, star_z], 
    1UnitfulAstro.Myr, 
    ["sim1" "sim2"],
)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "compare_CMDFPlot.png"))

figure = birthHistogramPlot(birth_pos, bins = 50)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "birthHistogramPlot.png"))

figure = timeSeriesPlot(time_series, mass_factor = 0, number_factor = 4)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "timeSeriesPlot.png"))

figure = scaleFactorSeriesPlot(time_series, mass_factor = 0, number_factor = 4)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "scaleFactorSeriesPlot.png"))

figure = redshiftSeriesPlot(time_series, mass_factor = 0, number_factor = 4)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "redshiftSeriesPlot.png"))

figure = compareSimulationsPlot(
    [time_series, time_series],
    "star_mass",
    "sfr",
    ["sim1" "sim2"];
    title = "SFMS relation",
    x_factor = 10,
    scale = [:identity, :log10],
)
Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "compareSimulationsPlot.png"),
)

figure = densityHistogramPlot(density, 1UnitfulAstro.Myr, factor = 10)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "densityHistogramPlot.png"))

figure = densityProfilePlot(
    pos, 
    gas_mass, 
    1UnitfulAstro.Myr, 
    scale = :log10, 
    bins = 50, 
    factor = 6,
)
Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "densityProfilePlot.png"),
)

figure = densityProfilePlot(
    [pos, pos],
    [gas_mass, gas_mass],
    1UnitfulAstro.Myr,
    ["sim_1" "sim_2"],
    bins = 50,
    scale = :log10,
    factor = 6,
)
Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "compare_densityProfilePlot.png"),
)

figure = metallicityProfilePlot(
    pos, 
    gas_mass, 
    gas_z, 
    1UnitfulAstro.Myr, 
    scale = :log10, 
    bins = 50,
)
Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "metallicityProfilePlot.png"),
)

figure = metallicityProfilePlot(
    [pos, pos],
    [gas_mass, gas_mass],
    [gas_z, gas_z],
    1UnitfulAstro.Myr,
    ["sim_1" "sim_2"],
    scale = :log10,
    bins = 50,
)
Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "compare_metallicityProfilePlot.png"),
)

figure = massProfilePlot(
    pos, 
    gas_mass, 
    1UnitfulAstro.Myr, 
    scale = :log10, 
    bins = 50, 
    factor = 10,
)
Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "massProfilePlot.png"),
)

figure = massProfilePlot(
    [pos, pos],
    [gas_mass, gas_mass],
    1UnitfulAstro.Myr,
    ["sim_1" "sim_2"],
    scale = :log10,
    bins = 50,
    factor = 10,
)
Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "compare_massProfilePlot.png"),
)

figure = sfrTxtPlot(
    sfrtxt_data,
    1,
    [4, 6],
    title = "run_A_01",
    bins = 50,
    scale = (:identity, :log10),
)
Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "compare_columns_sfrTxtPlot.png"),
)

figure = sfrTxtPlot(
    [sfrtxt_data, sfrtxt_data],
    1,
    6,
    ["sim_1" "sim_2"],
    title = "Column 6 vs 1",
    bins = 50,
    scale = (:identity, :log10),
)
Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "compare_sims_sfrTxtPlot.png"),
)

figure = temperatureHistogramPlot(temp_data, 1UnitfulAstro.Myr, bins = 30)
Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "temperatureHistogramPlot.png"),
)

figure = rhoTempPlot(temp_data, density, 1UnitfulAstro.Myr)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "rhoTempPlot.png"))

figure = KennicuttSchmidtPlot(
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
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "KennicuttSchmidtPlot.png"))