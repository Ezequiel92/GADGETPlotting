############################################################################################
# Testing script for GADGETPlotting.jl.
# It can be run as is and it shouldn't throw any errors.
#
# NOTE:
# This script is not intended to be a comprehensive test of the functionality, 
# but just to be a check that nothing is seriously broken, and to give examples
# on how to use the functions.
############################################################################################

include("src/GADGETPlotting.jl")

"Base path for the directories where the figures and animations will be saved."
const BASE_OUT_PATH = "results/"

"Directory containing the snapshot files."
const SNAP_PATH = "test_snapshots/"

"Base name of the snapshot files."
const SNAP_NAME = "snap"

"Side dimension of the simulated region, for the case of vacuum boundary conditions."
const BOX_SIZE = 200UnitfulAstro.kpc

"Value of ComovingIntegrationOn: 0 -> Newtonian simulation, 1 -> Cosmological simulation."
const SIM_COSMO = 0

"Frame rate for the animations."
const FPS = 4

"Particular snapshot index for testing purposes."
const SNAP_N = 21

mkpath(BASE_OUT_PATH)

############################################################################################
# TEST OF DATA ACQUISITION FUNCTIONS.
############################################################################################

snaps = getSnapshots(SNAP_NAME, SNAP_PATH)
snap_files = snaps["snap_files"]
# display(snaps)
# println()

# time_series = timeSeriesData(snap_files, sim_cosmo = SIM_COSMO)
# display(time_series)
# println()

pos = positionData(
    snap_files[SNAP_N],
    sim_cosmo = SIM_COSMO,
    box_size = BOX_SIZE,
    length_unit = UnitfulAstro.Mpc,
)
# display(pos)
# println()

# density = densityData(snap_files[SNAP_N], sim_cosmo = SIM_COSMO)
# display(density)
# println()

# hsml = hsmlData(snap_files[SNAP_N], sim_cosmo = SIM_COSMO)
# display(hsml)
# println()

gas_mass = massData(snap_files[SNAP_N], "gas", sim_cosmo = SIM_COSMO)
# display(gas_mass)
# println()

# dm_mass = massData(snap_files[SNAP_N], "dark_matter", sim_cosmo = SIM_COSMO)
# display(dm_mass)
# println()

star_mass = massData(snap_files[SNAP_N], "stars", sim_cosmo = SIM_COSMO)
# display(star_mass)
# println()

gas_z = zData(snap_files[SNAP_N], "gas", sim_cosmo = SIM_COSMO)
# display(gas_z)
# println()

star_z = zData(snap_files[SNAP_N], "stars", sim_cosmo = SIM_COSMO)
# display(star_z)
# println()

# sfrtxtdata = sfrTxtData(SNAP_PATH, SNAP_NAME, sim_cosmo = SIM_COSMO)
# display(sfrtxtdata)
# println()

# birth_pos = birthPlace(SNAP_N, snap_files, time_series["clock_time"], sim_cosmo = SIM_COSMO)
# display(birth_pos["birth_place"])
# println()

# ############################################################################################
# # TEST OF PLOTTING FUNCTIONS.
# ############################################################################################

# scatterGridPlot(pos)
# savefig(BASE_OUT_PATH * "test_scatterGridPlot.png")

# densityMapPlot(pos, gas_mass, density, hsml, plane = "XY")
# savefig(BASE_OUT_PATH * "test_densityMapPlot_XY.png")

# densityMapPlot(pos, gas_mass, density, hsml, axes = true)
# savefig(BASE_OUT_PATH * "test_densityMapPlot_All.png")

# color = [:batlow, :bone, :CMRmap, :grayC, :inferno, :seaborn_rocket_gradient, :YlOrRd_9]

# for c in color
#     densityMapPlot(pos, gas_mass, density, hsml, color = c, axes = true)
#     savefig(BASE_OUT_PATH * "test_densityMapPlot_" * string(c) * ".png")
# end

# starMapPlot(pos, plane = "XY", axes = true)
# savefig(BASE_OUT_PATH * "test_starMapPlot_XY.png")

# starMapPlot(pos, axes = true)
# savefig(BASE_OUT_PATH * "test_starMapPlot_All.png")

# gasStarEvolutionPlot(time_series, pos, SNAP_N)
# savefig(BASE_OUT_PATH * "test_gasStarEvolutionPlot.png")

# pgfplotsx()

# CMDFPlot(star_mass, star_z)
# savefig(BASE_OUT_PATH * "test_CMDFPlot.png")

# birthHistogramPlot(birth_pos, bins = 50)
# savefig(BASE_OUT_PATH * "test_birthHistogramPlot.png")

# timeSeriesPlot(time_series, mass_factor = 0, number_factor = 4)
# savefig(BASE_OUT_PATH * "test_timeSeriesPlot.png")

# scaleFactorSeriesPlot(time_series, mass_factor = 0, number_factor = 4)
# savefig(BASE_OUT_PATH * "test_scaleFactorSeriesPlot.png")

# redshiftSeriesPlot(time_series, mass_factor = 0, number_factor = 4)
# savefig(BASE_OUT_PATH * "test_redshiftSeriesPlot.png")

# compareSimulationsPlot(
#     [time_series, time_series],
#     "star_mass",
#     "sfr",
#     ["sim1" "sim2"];
#     title = "SFMS relation",
#     x_factor = 10,
#     scale = :log10,
# )
# savefig(BASE_OUT_PATH * "test_compareSimulationsPlot.png")

# densityHistogramPlot(density, 1UnitfulAstro.Myr, factor = 10)
# savefig(BASE_OUT_PATH * "test_densityHistogramPlot.png")

# densityProfilePlot(pos, star_mass, 1UnitfulAstro.Myr, scale = :log10, bins = 50)
# savefig(BASE_OUT_PATH * "test_star_densityProfilePlot.png")

# densityProfilePlot(pos, gas_mass, 1UnitfulAstro.Myr, scale = :log10, bins = 50, factor = 6)
# savefig(BASE_OUT_PATH * "test_gas_densityProfilePlot.png")

# densityProfilePlot(
#     [pos, pos],
#     [gas_mass, gas_mass],
#     1UnitfulAstro.Myr,
#     ["sim_1" "sim_2"],
#     bins = 50,
#     scale = :log10,
#     factor = 6,
# )
# savefig(BASE_OUT_PATH * "test_compare_gas_densityProfilePlot.png")

# densityProfilePlot(
#     [pos, pos],
#     [star_mass, star_mass],
#     1UnitfulAstro.Myr,
#     ["sim_1" "sim_2"],
#     bins = 50,
#     scale = :log10,
# )
# savefig(BASE_OUT_PATH * "test_compare_star_densityProfilePlot.png")

# metallicityProfilePlot(pos, gas_mass, gas_z, 1UnitfulAstro.Myr, scale = :log10, bins = 50)
# savefig(BASE_OUT_PATH * "test_gas_metallicityProfilePlot.png")

# metallicityProfilePlot(pos, star_mass, star_z, 1UnitfulAstro.Myr, scale = :log10, bins = 50)
# savefig(BASE_OUT_PATH * "test_star_metallicityProfilePlot.png")

# metallicityProfilePlot(
#     [pos, pos],
#     [gas_mass, gas_mass],
#     [gas_z, gas_z],
#     1UnitfulAstro.Myr,
#     ["sim_1" "sim_2"],
#     scale = :log10,
#     bins = 50,
# )
# savefig(BASE_OUT_PATH * "test_compare_gas_metallicityProfilePlot.png")

# metallicityProfilePlot(
#     [pos, pos],
#     [star_mass, star_mass],
#     [star_z, star_z],
#     1UnitfulAstro.Myr,
#     ["sim_1" "sim_2"],
#     scale = :log10,
#     bins = 50,
# )
# savefig(BASE_OUT_PATH * "test_compare_star_metallicityProfilePlot.png")

# massProfilePlot(pos, star_mass, 1UnitfulAstro.Myr, scale = :log10, bins = 50)
# savefig(BASE_OUT_PATH * "test_star_massProfilePlot.png")

# massProfilePlot(pos, gas_mass, 1UnitfulAstro.Myr, scale = :log10, bins = 50, factor = 10)
# savefig(BASE_OUT_PATH * "test_gas_massProfilePlot.png")

# massProfilePlot(
#     [pos, pos],
#     [gas_mass, gas_mass],
#     1UnitfulAstro.Myr,
#     ["sim_1" "sim_2"],
#     scale = :log10,
#     bins = 50,
#     factor = 10,
# )
# savefig(BASE_OUT_PATH * "test_compare_gas_massProfilePlot.png")

# massProfilePlot(
#     [pos, pos],
#     [star_mass, star_mass],
#     1UnitfulAstro.Myr,
#     ["sim_1" "sim_2"],
#     scale = :log10,
#     bins = 50,
# )
# savefig(BASE_OUT_PATH * "test_compare_star_massProfilePlot.png")

# sfr_data = sfrTxtData(SNAP_PATH, SNAP_NAME; sim_cosmo = SIM_COSMO)
# sfrTxtPlot(sfr_data, 1, [4, 6], title = "run_A_01", bins = 50, scale = (:identity, :log10))
# savefig(BASE_OUT_PATH * "test_sfrTxtPlot.png")

# ############################################################################################
# # TEST OF PIPELINE FUNCTIONS.
# ############################################################################################

# gr()

# scatterGridPipeline(
#     SNAP_NAME,
#     SNAP_PATH,
#     "scatter_animation",
#     FPS,
#     output_path = BASE_OUT_PATH * "scatter_grid/",
#     sim_cosmo = SIM_COSMO,
#     box_size = BOX_SIZE,
#     length_unit = UnitfulAstro.Mpc,
# )

# densityMapPipeline(
#     SNAP_NAME,
#     SNAP_PATH,
#     "density_animation",
#     FPS,
#     output_path = BASE_OUT_PATH * "density_map/",
#     sim_cosmo = SIM_COSMO,
#     plane = "All",
#     box_size = BOX_SIZE,
# )

# starMapPipeline(
#     SNAP_NAME,
#     SNAP_PATH,
#     "star_animation",
#     FPS,
#     output_path = BASE_OUT_PATH * "star_map/",
#     sim_cosmo = SIM_COSMO,
#     plane = "All",
#     box_size = BOX_SIZE,
# )

# gasStarEvolutionPipeline(
#     SNAP_NAME,
#     SNAP_PATH,
#     "gas_star_evolution",
#     FPS,
#     output_path = BASE_OUT_PATH * "gas_star_evolution/",
#     sim_cosmo = SIM_COSMO,
#     box_size = BOX_SIZE,
# )

pgfplotsx()

# evolutionSummaryPipeline(
#     SNAP_NAME,
#     SNAP_PATH,
#     "evolution_summary",
#     output_path = BASE_OUT_PATH * "evolution_summary/",
#     sim_cosmo = SIM_COSMO,
#     mass_factor = 10,
#     number_factor = 4,
# )

# compareSimulationsPipeline(
#     [SNAP_NAME, SNAP_NAME],
#     [SNAP_PATH, SNAP_PATH],
#     ["sim1" "sim2"],
#     "compare",
#     "star_mass",
#     "sfr",
#     output_path = BASE_OUT_PATH * "compare_simulations/",
#     sim_cosmo = SIM_COSMO,
#     title = "SFMS relation",
#     x_factor = 10,
# )

# densityHistogramPipeline(
#     SNAP_NAME,
#     SNAP_PATH,
#     "density_histogram_animation",
#     FPS,
#     output_path = BASE_OUT_PATH * "density_histogram/",
#     sim_cosmo = SIM_COSMO,
#     factor = 10,
# )

# densityProfilePipeline(
#     SNAP_NAME,
#     SNAP_PATH,
#     "density_profile_animation",
#     FPS,
#     "stars";
#     output_path = BASE_OUT_PATH * "stars_density_profile/",
#     sim_cosmo = SIM_COSMO,
#     scale = :log10,
#     bins = 80,
#     factor = 6,
#     box_factor = 0.125,
#     box_size = BOX_SIZE,
# )

# densityProfilePipeline(
#     [SNAP_NAME, SNAP_NAME],
#     [SNAP_PATH, SNAP_PATH],
#     "compare_density_profile_animation",
#     FPS,
#     "gas",
#     ["sim1" "sim2"];
#     output_path = BASE_OUT_PATH * "compare_gas_density_profile/",
#     sim_cosmo = SIM_COSMO,
#     scale = :log10,
#     bins = 60,
#     factor = 6,
#     box_factor = 0.125,
#     box_size = BOX_SIZE,
# )

# metallicityProfilePipeline(
#     SNAP_NAME,
#     SNAP_PATH,
#     "metallicity_profile_animation",
#     FPS,
#     "stars";
#     output_path = BASE_OUT_PATH * "stars_metallicity_profile/",
#     sim_cosmo = SIM_COSMO,
#     scale = :log10,
#     bins = 80,
#     box_factor = 5.0,
#     box_size = BOX_SIZE,
# )

# metallicityProfilePipeline(
#     [SNAP_NAME, SNAP_NAME],
#     [SNAP_PATH, SNAP_PATH],
#     "compare_metallicity_profile_animation",
#     FPS,
#     "gas",
#     ["sim1" "sim2"];
#     output_path = BASE_OUT_PATH * "compare_gas_metallicity_profile/",
#     sim_cosmo = SIM_COSMO,
#     scale = :log10,
#     bins = 60,
#     box_factor = 5.0,
#     box_size = BOX_SIZE,
# )

# massProfilePipeline(
#     SNAP_NAME,
#     SNAP_PATH,
#     "mass_profile_animation",
#     FPS,
#     "stars";
#     output_path = BASE_OUT_PATH * "stars_mass_profile/",
#     sim_cosmo = SIM_COSMO,
#     scale = :log10,
#     bins = 80,
#     factor = 6,
#     box_factor = 0.125,
#     box_size = BOX_SIZE,
# )

# massProfilePipeline(
#     [SNAP_NAME, SNAP_NAME],
#     [SNAP_PATH, SNAP_PATH],
#     "compare_mass_profile_animation",
#     FPS,
#     "gas",
#     ["sim1" "sim2"];
#     output_path = BASE_OUT_PATH * "compare_gas_mass_profile/",
#     sim_cosmo = SIM_COSMO,
#     scale = :log10,
#     bins = 60,
#     factor = 6,
#     box_factor = 0.125,
#     box_size = BOX_SIZE,
# )

# CMDFPipeline(
#     SNAP_NAME,
#     SNAP_PATH,
#     "CMDF_animation",
#     FPS,
#     output_path = BASE_OUT_PATH * "CMDF/",
#     sim_cosmo = SIM_COSMO,
# )

# birthHistogramPipeline(
#     SNAP_NAME,
#     SNAP_PATH,
#     "birth_histogram_animation",
#     FPS,
#     output_path = BASE_OUT_PATH * "birth_histogram/",
#     sim_cosmo = SIM_COSMO,
# )

sfrTxtPipeline(
    [SNAP_NAME, SNAP_NAME], 
    [SNAP_PATH, SNAP_PATH], 
    1, 
    [4, 6], 
    output_path = BASE_OUT_PATH * "sfr_txt/",
    sim_cosmo = SIM_COSMO, 
    title = ["sim_1", "sim_2"], 
    bins = 50, 
    scale = (:identity, :log10),
)

############################################################################################
# TEST OF AUXILIARY FUNCTIONS.
############################################################################################

fig2D = plot(rand(100))
println(relative(fig2D, 0.5, 0.5))
fig3D = surface(rand(100, 100))
println(relative(fig3D, 0.5, 0.5, 0.5))

makeVideo(BASE_OUT_PATH * "scatter_grid/images", ".png", BASE_OUT_PATH, "test_video", FPS)

x_data = [1:1000...]
y_data = rand(1000)
x_smooth, y_smooth = smoothWindow(x_data, y_data, 50)
plot(x_data, y_data, seriestype = :scatter, legend = false)
plot!(x_smooth, y_smooth, seriestype = :line, lw = 2)
savefig(BASE_OUT_PATH * "test_smoothWindow.png")

positions = pos["gas"]
distances = sqrt.(positions[1, :] .^ 2 .+ positions[2, :] .^ 2 .+ positions[3, :] .^ 2)
box_size = ustrip(Float64, UnitfulAstro.Mpc, BOX_SIZE)

r, ρ = densityProfile(gas_mass["mass"], distances, box_size, 80)
plot(r, ρ, lw = 2, xlabel = "r / $(pos["unit"])", ylabel = L"\rho", legend = false)
savefig(BASE_OUT_PATH * "test_densityProfile.png")

r, z = metallicityProfile(gas_mass["mass"], distances, gas_z["Z"], box_size, 80)
plot(r, z, lw = 2, xlabel = "r / $(pos["unit"])", ylabel = "Z / Zsun", legend = false)
savefig(BASE_OUT_PATH * "test_metallicityProfile.png")

r, m = massProfile(gas_mass["mass"], distances, box_size, 80)
plot(r, m, lw = 2, xlabel = "r / $(pos["unit"])", ylabel = "Mass", legend = false)
savefig(BASE_OUT_PATH * "test_massProfile.png")

max_z = findmax(star_z["Z"])
max_Z = max_z[1] / star_mass["mass"][max_z[2]]
z, m = CMDF(star_mass["mass"], star_z["Z"], max_Z, 50)
plot(z, m, 
    lw = 2, 
    xlabel = "Z", 
    ylabel = L"M_{\star}(< Z) \, / \, M_{\star}", 
    legend = false
    )
savefig(BASE_OUT_PATH * "test_CMDF.png")

println("Everything worked just fine!!")

############################################################################################
# DELETE ALL GENERATED TESTING FILES.
############################################################################################

# rm(BASE_OUT_PATH, recursive = true)