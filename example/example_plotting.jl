############################################################################################
# Example plotting functions
############################################################################################

figure = scatter_grid_plot(pos)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "scatter_grid_plot.svg"))

figure = density_map_plot(pos, gas_mass, density, hsml)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "density_map_plot.png"))

figure = star_map_plot(pos)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "star_map_plot.png"))

figure = gas_star_evolution_plot(SNAP_N, time_series, pos)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "gas_star_evolution_plot.png"))

figure = cmdf_plot(star_mass, star_z, 1UnitfulAstro.Myr)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "cmdf_plot.png"))

figure = cmdf_plot(
    [star_mass, star_mass], 
    [star_z, star_z], 
    1UnitfulAstro.Myr, 
    ["sim1" "sim2"],
)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "compare_cmdf_plot.png"))

figure = birth_histogram_plot(birth_pos, bins = 50)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "birth_histogram_plot.png"))

figure = time_series_plot(time_series, mass_factor = 0, number_factor = 4)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "time_series_plot.png"))

figure = scale_factor_series_plot(time_series, mass_factor = 0, number_factor = 4)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "scale_factor_series_plot.png"))

figure = redshift_series_plot(time_series, mass_factor = 0, number_factor = 4)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "redshift_series_plot.png"))

figure = compare_simulations_plot(
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
    joinpath(BASE_OUT_PATH, "compare_simulations_plot.png"),
)

figure = density_histogram_plot(density, 1UnitfulAstro.Myr, factor = 10)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "density_histogram_plot.png"))

figure = density_profile_plot(
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
    joinpath(BASE_OUT_PATH, "density_profile_plot.png"),
)

figure = density_profile_plot(
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
    joinpath(BASE_OUT_PATH, "compare_density_profile_plot.png"),
)

figure = metallicity_profile_plot(
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
    joinpath(BASE_OUT_PATH, "metallicity_profile_plot.png"),
)

figure = metallicity_profile_plot(
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
    joinpath(BASE_OUT_PATH, "compare_metallicity_profile_plot.png"),
)

figure = mass_profile_plot(
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
    joinpath(BASE_OUT_PATH, "mass_profile_plot.png"),
)

figure = mass_profile_plot(
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
    joinpath(BASE_OUT_PATH, "compare_mass_profile_plot.png"),
)

figure = sfr_txt_plot(
    sfr_txt_data,
    1,
    [4, 6],
    title = "run_A_01",
    bins = 50,
    scale = (:identity, :log10),
)
Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "compare_columns_sfr_txt_plot.png"),
)

figure = sfr_txt_plot(
    [sfr_txt_data, sfr_txt_data],
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
    joinpath(BASE_OUT_PATH, "compare_sims_sfr_txt_plot.png"),
)

figure = temperature_histogram_plot(temp_data, 1UnitfulAstro.Myr, bins = 30)
Base.invokelatest(
    savefig, 
    figure, 
    joinpath(BASE_OUT_PATH, "temperature_histogram_plot.png"),
)

figure = rho_temp_plot(temp_data, density, 1UnitfulAstro.Myr)
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "rho_temp_plot.png"))

figure = kennicutt_schmidt_plot(
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
Base.invokelatest(savefig, figure, joinpath(BASE_OUT_PATH, "kennicutt_schmidt_plot.png"))