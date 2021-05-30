############################################################################################
# Example pipeline functions
############################################################################################

scatter_grid_pipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "scatter_animation",
    FPS,
    output_path = joinpath(BASE_OUT_PATH, "scatter_grid"),
    sim_cosmo = SIM_COSMO,
    box_size = BOX_SIZE,
    length_unit = UnitfulAstro.Mpc,
)

density_map_pipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "density_animation",
    FPS,
    output_path = joinpath(BASE_OUT_PATH, "density_map"),
    sim_cosmo = SIM_COSMO,
    plane = "All",
    box_size = BOX_SIZE,
)

star_map_pipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "star_animation",
    FPS,
    output_path = joinpath(BASE_OUT_PATH, "star_map"),
    sim_cosmo = SIM_COSMO,
    plane = "All",
    box_size = BOX_SIZE,
)

gas_star_evolution_pipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "gas_star_evolution",
    FPS,
    output_path = joinpath(BASE_OUT_PATH, "gas_star_evolution"),
    sim_cosmo = SIM_COSMO,
    box_size = BOX_SIZE,
)

cmdf_pipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "CMDF_animation",
    FPS,
    output_path = joinpath(BASE_OUT_PATH, "CMDF"),
    sim_cosmo = SIM_COSMO,
)

cmdf_pipeline(
    [SNAP_NAME, SNAP_NAME],
    [BASE_SRC_PATH, BASE_SRC_PATH],
    "CMDF_animation",
    FPS,
    ["sim1" "sim2"],
    output_path = joinpath(BASE_OUT_PATH, "compare_CMDF"),
    sim_cosmo = SIM_COSMO,
)

birth_histogram_pipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "birth_histogram_animation",
    FPS,
    output_path = joinpath(BASE_OUT_PATH, "birth_histogram"),
    sim_cosmo = SIM_COSMO,
)

evolution_summary_pipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "evolution_summary",
    output_path = joinpath(BASE_OUT_PATH, "evolution_summary"),
    sim_cosmo = SIM_COSMO,
    mass_factor = 10,
    number_factor = 4,
)

compare_simulations_pipeline(
    [SNAP_NAME, SNAP_NAME],
    [BASE_SRC_PATH, BASE_SRC_PATH],
    ["sim1" "sim2"],
    "compare",
    "star_mass",
    "sfr",
    output_path = joinpath(BASE_OUT_PATH, "compare_simulations"),
    sim_cosmo = SIM_COSMO,
    title = "SFMS relation",
    x_factor = 10,
    text_quantity = "star_mass"
)

density_histogram_pipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "density_histogram_animation",
    FPS,
    output_path = joinpath(BASE_OUT_PATH, "density_histogram"),
    sim_cosmo = SIM_COSMO,
    factor = 10,
)

density_profile_pipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "density_profile_animation",
    FPS,
    "stars";
    output_path = joinpath(BASE_OUT_PATH, "stars_density_profile"),
    sim_cosmo = SIM_COSMO,
    scale = :log10,
    bins = 80,
    factor = 6,
    box_factor = 0.125,
    box_size = BOX_SIZE,
)

density_profile_pipeline(
    [SNAP_NAME, SNAP_NAME],
    [BASE_SRC_PATH, BASE_SRC_PATH],
    "compare_density_profile_animation",
    FPS,
    "gas",
    ["sim1" "sim2"];
    output_path = joinpath(BASE_OUT_PATH, "compare_gas_density_profile"),
    sim_cosmo = SIM_COSMO,
    scale = :log10,
    bins = 60,
    factor = 6,
    box_factor = 0.125,
    box_size = BOX_SIZE,
)

metallicity_profile_pipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "metallicity_profile_animation",
    FPS,
    "stars";
    output_path = joinpath(BASE_OUT_PATH, "stars_metallicity_profile"),
    sim_cosmo = SIM_COSMO,
    scale = :log10,
    bins = 80,
    box_factor = 5.0,
    box_size = BOX_SIZE,
)

metallicity_profile_pipeline(
    [SNAP_NAME, SNAP_NAME],
    [BASE_SRC_PATH, BASE_SRC_PATH],
    "compare_metallicity_profile_animation",
    FPS,
    "gas",
    ["sim1" "sim2"];
    output_path = joinpath(BASE_OUT_PATH, "compare_gas_metallicity_profile"),
    sim_cosmo = SIM_COSMO,
    scale = :log10,
    bins = 60,
    box_factor = 5.0,
    box_size = BOX_SIZE,
)

mass_profile_pipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "mass_profile_animation",
    FPS,
    "stars";
    output_path = joinpath(BASE_OUT_PATH, "stars_mass_profile"),
    sim_cosmo = SIM_COSMO,
    scale = :log10,
    bins = 80,
    factor = 6,
    box_factor = 0.125,
    box_size = BOX_SIZE,
)

mass_profile_pipeline(
    [SNAP_NAME, SNAP_NAME],
    [BASE_SRC_PATH, BASE_SRC_PATH],
    "compare_mass_profile_animation",
    FPS,
    "gas",
    ["sim1" "sim2"];
    output_path = joinpath(BASE_OUT_PATH, "compare_gas_mass_profile"),
    sim_cosmo = SIM_COSMO,
    scale = :log10,
    bins = 60,
    factor = 6,
    box_factor = 0.125,
    box_size = BOX_SIZE,
)

sfr_txt_pipeline(
    [SNAP_NAME * "_000", SNAP_NAME * "_000"],
    [BASE_SRC_PATH, BASE_SRC_PATH],
    1,
    [4, 6],
    output_path = joinpath(BASE_OUT_PATH, "sfr_txt_column_compare"),
    sim_cosmo = SIM_COSMO,
    titles = ["sim_1", "sim_2"],
    bins = 50,
    scale = (:identity, :log10),
)

sfr_txt_pipeline(
    [SNAP_NAME * "_000", SNAP_NAME * "_000"],
    [BASE_SRC_PATH, BASE_SRC_PATH],
    1,
    [4, 6],
    output_path = joinpath(BASE_OUT_PATH, "sfr_txt_sim_compare"),
    sim_cosmo = SIM_COSMO,
    comparison_type = 1,
    titles = ["Column 4 vs. 1", "Column 6 vs. 1"],
    labels = ["sim_1" "sim_2"],
    bins = 50,
    scale = (:identity, :log10),
)

temperature_histogram_pipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "temperature_histogram_animation",
    FPS,
    output_path = joinpath(BASE_OUT_PATH, "temperature_histogram"),
    sim_cosmo = SIM_COSMO,
)

rho_temp_pipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "rho_vs_temp_animation",
    FPS,
    output_path = joinpath(BASE_OUT_PATH, "rho_vs_temp"),
    sim_cosmo = SIM_COSMO,
)

kennicutt_schmidt_pipeline(
    SNAP_NAME,
    BASE_SRC_PATH;
    output_path = joinpath(BASE_OUT_PATH, "Kennicutt_Schmidt"),
    sim_cosmo = SIM_COSMO,
    max_r = BOX_SIZE,
    bins = 80,
    error_formating = "conf_interval",
    time_unit = UnitfulAstro.yr,
)

# The cpu.txt of cosmological simulations is too heavy for GitHub, and `get_cpu_txt`
# doesn't distinguish among different types of simulations.
if SIM_COSMO == 0
    cpu_txt_pipeline(
        [BASE_SRC_PATH, BASE_SRC_PATH],
        ["i/o", "hotngbs", "density"];
        output_path = joinpath(BASE_OUT_PATH, "cpu_txt"),
        titles = ["title_1",],
        names = ["sim_1",],
    )

    cpu_txt_pipeline(
        [BASE_SRC_PATH, BASE_SRC_PATH],
        "density",
        ["sim_1" "sim_2"];
        output_path = joinpath(BASE_OUT_PATH, "compare_cpu_txt"),
        title = "title",
    )
end