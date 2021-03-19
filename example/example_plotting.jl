############################################################################################
# PLOTTING FUNCTIONS.
############################################################################################

figure = scatterGridPlot(pos)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "test_scatterGridPlot.png"))

figure = densityMapPlot(pos, gas_mass, density, hsml, plane = "XY")
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "test_densityMapPlot_XY.png"))

figure = densityMapPlot(pos, gas_mass, density, hsml, axes = true)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "test_densityMapPlot_All.png"))

color = [:batlow, :bone, :CMRmap, :grayC, :inferno, :seaborn_rocket_gradient, :YlOrRd_9]

for c in color
    figure = densityMapPlot(pos, gas_mass, density, hsml, color = c, axes = true)
    Base.invokelatest(
        savefig, 
        figure, 
        joinpath(BASE_OUT_PATH, "test_densityMapPlot_" * string(c) * ".png"),
    )
end

figure = starMapPlot(pos, plane = "XY", axes = true)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "test_starMapPlot_XY.png"))

figure = starMapPlot(pos, axes = true)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "test_starMapPlot_All.png"))

figure = gasStarEvolutionPlot(SNAP_N, time_series, pos)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "test_gasStarEvolutionPlot.png"))

figure = CMDFPlot(star_mass, star_z, 1UnitfulAstro.Myr)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "test_CMDFPlot.png"))

figure = CMDFPlot(
    [star_mass, star_mass], 
    [star_z, star_z], 
    1UnitfulAstro.Myr, 
    ["sim1" "sim2"],
)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "test_compare_CMDFPlot.png"))

figure = birthHistogramPlot(birth_pos, bins = 50)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "test_birthHistogramPlot.png"))

figure = timeSeriesPlot(time_series, mass_factor = 0, number_factor = 4)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "test_timeSeriesPlot.png"))

figure = scaleFactorSeriesPlot(time_series, mass_factor = 0, number_factor = 4)
Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "test_scaleFactorSeriesPlot.png"),
)

figure = redshiftSeriesPlot(time_series, mass_factor = 0, number_factor = 4)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "test_redshiftSeriesPlot.png"))

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
    joinpath(BASE_OUT_PATH, "test_compareSimulationsPlot.png"),
)

figure = densityHistogramPlot(density, 1UnitfulAstro.Myr, factor = 10)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "test_densityHistogramPlot.png"))

figure = densityProfilePlot(pos, star_mass, 1UnitfulAstro.Myr, scale = :log10, bins = 50)
Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "test_star_densityProfilePlot.png"),
)

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
    joinpath(BASE_OUT_PATH, "test_gas_densityProfilePlot.png"),
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
    joinpath(BASE_OUT_PATH, "test_compare_gas_densityProfilePlot.png"),
)

figure = densityProfilePlot(
    [pos, pos],
    [star_mass, star_mass],
    1UnitfulAstro.Myr,
    ["sim_1" "sim_2"],
    bins = 50,
    scale = :log10,
)
Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "test_compare_star_densityProfilePlot.png"),
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
    joinpath(BASE_OUT_PATH, "test_gas_metallicityProfilePlot.png"),
)

figure = metallicityProfilePlot(
    pos, 
    star_mass, 
    star_z, 
    1UnitfulAstro.Myr, 
    scale = :log10, 
    bins = 50,
)
Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "test_star_metallicityProfilePlot.png"),
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
    joinpath(BASE_OUT_PATH, "test_compare_gas_metallicityProfilePlot.png"),
)

figure = metallicityProfilePlot(
    [pos, pos],
    [star_mass, star_mass],
    [star_z, star_z],
    1UnitfulAstro.Myr,
    ["sim_1" "sim_2"],
    scale = :log10,
    bins = 50,
)
figure = Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "test_compare_star_metallicityProfilePlot.png"),
)

figure = massProfilePlot(pos, star_mass, 1UnitfulAstro.Myr, scale = :log10, bins = 50)
Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "test_star_massProfilePlot.png"),
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
    joinpath(BASE_OUT_PATH, "test_gas_massProfilePlot.png"),
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
    joinpath(BASE_OUT_PATH, "test_compare_gas_massProfilePlot.png"),
)

figure = massProfilePlot(
    [pos, pos],
    [star_mass, star_mass],
    1UnitfulAstro.Myr,
    ["sim_1" "sim_2"],
    scale = :log10,
    bins = 50,
)
Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "test_compare_star_massProfilePlot.png"),
)

figure = sfrTxtPlot(
    sfrtxt_data,
    1,
    [4, 6],
    title = "run_A_01",
    bins = 50,
    scale = (:identity, :log10),
)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "test_sfrTxtPlot.png"))

figure = temperatureHistogramPlot(temp_data, 1UnitfulAstro.Myr, bins = 30)
Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "test_temperatureHistogramPlot.png"),
)

figure = rhoTempPlot(temp_data, density, 1UnitfulAstro.Myr)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "test_rhoTempPlot.png"))

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
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "test_KennicuttSchmidtPlot.png"))