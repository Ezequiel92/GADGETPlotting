############################################################################################
# Test plotting functions
############################################################################################

@testset "Plotting functions" begin

    temp_img = joinpath(@__DIR__, "test_img.png")

    fig = @test_nowarn scatter_grid_plot(pos)
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "scatter_grid_plot.png") load(temp_img)

    # Commented out because the @info in sphMapping() is interpreted as an error output.
    # Until Ludwing corrects this detail this will not pass the tests.
    # fig = @test_nowarn density_map_plot(nothing, pos, gas_mass, density, hsml)
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "density_map_plot.png") load(temp_img)

    fig = @test_nowarn star_map_plot(pos)
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "star_map_plot_All.png") load(temp_img)

    fig = @test_nowarn gas_star_evolution_plot(SNAP_N, time_series, pos)
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "gas_star_evolution_plot.png") load(temp_img)

    fig = @test_nowarn cmdf_plot(star_mass, star_z, 1UnitfulAstro.Myr)
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "cmdf_plot.png") load(temp_img)

    fig = @test_nowarn cmdf_plot(
        [star_mass, star_mass], 
        [star_z, star_z], 
        1UnitfulAstro.Myr, 
        ["sim1" "sim2"],
    )
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "compare_cmdf_plot.png") load(temp_img)

    fig = @test_nowarn birth_histogram_plot(birth_pos, bins = 50)
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "birth_histogram_plot.png") load(temp_img)

    fig = @test_nowarn time_series_plot(time_series, mass_factor = 0, number_factor = 4)
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "time_series_plot.png") load(temp_img)

    fig = @test_nowarn scale_factor_series_plot(
        time_series, 
        mass_factor = 0, 
        number_factor = 4,
    )
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "scale_factor_series_plot.png") load(temp_img)

    fig = @test_nowarn redshift_series_plot(time_series, mass_factor = 0, number_factor = 4)
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "redshift_series_plot.png") load(temp_img)

    fig = @test_nowarn compare_simulations_plot(
        [time_series, time_series],
        "star_mass",
        "sfr",
        ["sim1" "sim2"];
        title = "SFMS relation",
        x_factor = 10,
        scale = (:identity, :log10),
    )
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "compare_simulations_plot.png") load(temp_img)

    fig = @test_nowarn density_histogram_plot(density, 1UnitfulAstro.Myr, factor = 10)
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "density_histogram_plot.png") load(temp_img)

    fig = @test_nowarn fraction_histogram_plot(
        fmol, 
        1UnitfulAstro.Myr, 
        "Molecular fraction", 
        y_scale = :log10,
    )
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "fmol_histogram_plot.png") load(temp_img)

    fig = @test_nowarn fraction_histogram_plot(
        fatom, 
        1UnitfulAstro.Myr, 
        "Atomic fraction", 
        y_scale = :log10,
    )
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "fatom_histogram_plot.png") load(temp_img)

    fig = @test_nowarn density_profile_plot(
        pos, 
        gas_mass, 
        1UnitfulAstro.Myr, 
        scale = :log10, 
        bins = 50, 
        factor = 6,
    )
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "density_profile_plot.png") load(temp_img)
    
    fig = @test_nowarn density_profile_plot(
        [pos, pos],
        [gas_mass, gas_mass],
        1UnitfulAstro.Myr,
        ["sim_1" "sim_2"],
        bins = 50,
        scale = :log10,
        factor = 6,
    )
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "compare_density_profile_plot.png") load(temp_img)

    fig = @test_nowarn metallicity_profile_plot(
        pos, 
        gas_mass, 
        gas_z, 
        1UnitfulAstro.Myr, 
        scale = :log10, 
        bins = 50,
    )
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "metallicity_profile_plot.png") load(temp_img)

    fig = @test_nowarn metallicity_profile_plot(
        [pos, pos],
        [gas_mass, gas_mass],
        [gas_z, gas_z],
        1UnitfulAstro.Myr,
        ["sim_1" "sim_2"],
        scale = :log10,
        bins = 50,
    )
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "compare_metallicity_profile_plot.png") load(temp_img)

    fig = @test_nowarn mass_profile_plot(
        pos, 
        gas_mass, 
        1UnitfulAstro.Myr, 
        scale = :log10, 
        bins = 50, 
        factor = 10,
    )
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "mass_profile_plot.png") load(temp_img)

    fig = @test_nowarn mass_profile_plot(
        [pos, pos],
        [gas_mass, gas_mass],
        1UnitfulAstro.Myr,
        ["sim_1" "sim_2"],
        scale = :log10,
        bins = 50,
        factor = 10,
    )
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "compare_mass_profile_plot.png") load(temp_img)

    fig = @test_nowarn sfr_txt_plot(
        sfr_txt_data,
        1,
        [4, 6],
        title = "run_A_01",
        bins = 50,
        scale = (:identity, :log10),
    )
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "compare_columns_sfr_txt_plot.png") load(temp_img)

    fig = @test_nowarn sfr_txt_plot(
        [sfr_txt_data, sfr_txt_data],
        1,
        6,
        ["sim_1" "sim_2"],
        title = "Column 6 vs 1",
        bins = 50,
        scale = (:identity, :log10),
    )
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "compare_sims_sfr_txt_plot.png") load(temp_img)

    fig = @test_nowarn temperature_histogram_plot(temp_data, 1UnitfulAstro.Myr, bins = 30)
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "temperature_histogram_plot.png") load(temp_img)

    fig = @test_nowarn rho_temp_plot(temp_data, density, 1UnitfulAstro.Myr)
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "rho_temp_plot.png") load(temp_img)

    fig = @test_nowarn fraction_temp_plot(temp_data, fmol, 1UnitfulAstro.Myr, "Molecular fraction")
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "mol_fraction_temp_plot.png") load(temp_img)

    fig = @test_nowarn fraction_temp_plot(temp_data, fatom, 1UnitfulAstro.Myr, "Atomic fraction")
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "atom_fraction_temp_plot.png") load(temp_img)

    fig = @test_nowarn kennicutt_schmidt_plot(
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
    # savefig(fig, temp_img)
    # @test_reference joinpath(BASE_DATA_PATH, "kennicutt_schmidt_plot.png") load(temp_img)

    # The cpu.txt of cosmological simulations is too heavy for GitHub, and `get_cpu_txt`
    # doesn't distinguish among different types of simulations.
    if SIM_COSMO == 0
        figure = @test_nowarn cpu_txt_plot(cpu_txt_data)
        # savefig(fig, temp_img)
        # @test_reference joinpath(BASE_DATA_PATH, "cpu_txt_plot.png") load(temp_img)

        figure = @test_nowarn cpu_txt_plot([cpu_txt_data, cpu_txt_data], ["sim_1" "sim_2"])
        # savefig(fig, temp_img)
        # @test_reference joinpath(BASE_DATA_PATH, "compare_cpu_txt_plot.png") load(temp_img)
    end

    # `FMOL` is not a block present in the examples of cosmological snapshots
    if sim_cosmo == 0
        figure = @test_nowarn quantities_2D_plot(
            quantities2D,
            "SFE",
            "OH",
            Dict(
                "mass" => UnitfulAstro.Msun, 
                "length" => UnitfulAstro.Mpc, 
                "time" => UnitfulAstro.Myr,
            );
            title = "OH vs. SFE",
            x_factor = 28,
            y_factor = 0,
            scale = (:log10, :log10),
        )
        # savefig(fig, temp_img)
        # @test_reference joinpath(BASE_DATA_PATH, "quantities_2D_plot.png") load(temp_img)

        figure = @test_nowarn fatom_rho_plot(fatom, density, 1UnitfulAstro.Myr)
        # savefig(fig, temp_img)
        # @test_reference joinpath(BASE_DATA_PATH, "fatom_rho_plot.png") load(temp_img)

        figure = @test_nowarn fmol_Z_plot(fmol, gas_z, gas_mass, 1UnitfulAstro.Myr)
        # savefig(fig, temp_img)
        # @test_reference joinpath(BASE_DATA_PATH, "fmol_Z_plot.png") load(temp_img)

        figure = @test_nowarn fmol_fatom_plot(
            fmol, 
            fatom, 
            density, 
            gas_z, 
            gas_mass, 
            1UnitfulAstro.Myr,
        )
        # savefig(fig, temp_img)
        # @test_reference joinpath(BASE_DATA_PATH, "fmol_fatom_plot.png") load(temp_img)
    end
    
    # @test_nowarn rm(temp_img)

end