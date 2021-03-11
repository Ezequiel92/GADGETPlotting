############################################################################################
# Example script for GADGETPlotting.jl.
# When run as is, it shouldn't throw any errors.
#
# NOTE:
# This script is not intended to be a comprehensive test of the functionality, 
# but just to be a check that nothing is seriously broken, and to give an example
# on how to load GADGETPlotting.jl and how to use every functions.
############################################################################################

push!(LOAD_PATH, "./src/")
using GADGETPlotting, Unitful, UnitfulAstro, Plots, LaTeXStrings, GLM

"Base path to the directories where the output images and animations will be saved."
const BASE_OUT_PATH = joinpath(@__DIR__, "example_results/")

"Directory containing the snapshot files."
const BASE_SRC_PATH = joinpath(@__DIR__, "test_data/")

"Base name of the snapshot files, set in the GADGET variable SnapshotFileBase."
const SNAP_NAME = "snap"

"Side dimension of the simulated region, for the case of vacuum boundary conditions."
const BOX_SIZE = 200UnitfulAstro.kpc

"Value of ComovingIntegrationOn: 0 -> Newtonian simulation, 1 -> Cosmological simulation."
const SIM_COSMO = 0

"Frame rate for the animations."
const FPS = 4

"One particular snapshot index for testing purposes."
const SNAP_N = 21

mkpath(BASE_OUT_PATH)

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

############################################################################################
# PLOTTING FUNCTIONS.
############################################################################################

scatterGridPlot(pos)
savefig(BASE_OUT_PATH * "test_scatterGridPlot.png")

densityMapPlot(pos, gas_mass, density, hsml, plane = "XY")
savefig(BASE_OUT_PATH * "test_densityMapPlot_XY.png")

densityMapPlot(pos, gas_mass, density, hsml, axes = true)
savefig(BASE_OUT_PATH * "test_densityMapPlot_All.png")

color = [:batlow, :bone, :CMRmap, :grayC, :inferno, :seaborn_rocket_gradient, :YlOrRd_9]

for c in color
    densityMapPlot(pos, gas_mass, density, hsml, color = c, axes = true)
    savefig(BASE_OUT_PATH * "test_densityMapPlot_" * string(c) * ".png")
end

starMapPlot(pos, plane = "XY", axes = true)
savefig(BASE_OUT_PATH * "test_starMapPlot_XY.png")

starMapPlot(pos, axes = true)
savefig(BASE_OUT_PATH * "test_starMapPlot_All.png")

gasStarEvolutionPlot(time_series, pos, SNAP_N)
savefig(BASE_OUT_PATH * "test_gasStarEvolutionPlot.png")

pgfplotsx()

CMDFPlot(star_mass, star_z, 1UnitfulAstro.Myr)
savefig(BASE_OUT_PATH * "test_CMDFPlot.png")

CMDFPlot([star_mass, star_mass], [star_z, star_z], 1UnitfulAstro.Myr, ["sim1" "sim2"])
savefig(BASE_OUT_PATH * "test_compare_CMDFPlot.png")

birthHistogramPlot(birth_pos, bins = 50)
savefig(BASE_OUT_PATH * "test_birthHistogramPlot.png")

timeSeriesPlot(time_series, mass_factor = 0, number_factor = 4)
savefig(BASE_OUT_PATH * "test_timeSeriesPlot.png")

scaleFactorSeriesPlot(time_series, mass_factor = 0, number_factor = 4)
savefig(BASE_OUT_PATH * "test_scaleFactorSeriesPlot.png")

redshiftSeriesPlot(time_series, mass_factor = 0, number_factor = 4)
savefig(BASE_OUT_PATH * "test_redshiftSeriesPlot.png")

compareSimulationsPlot(
    [time_series, time_series],
    "star_mass",
    "sfr",
    ["sim1" "sim2"];
    title = "SFMS relation",
    x_factor = 10,
    scale = (:identity, :log10),
)
savefig(BASE_OUT_PATH * "test_compareSimulationsPlot.png")

densityHistogramPlot(density, 1UnitfulAstro.Myr, factor = 10)
savefig(BASE_OUT_PATH * "test_densityHistogramPlot.png")

densityProfilePlot(pos, star_mass, 1UnitfulAstro.Myr, scale = :log10, bins = 50)
savefig(BASE_OUT_PATH * "test_star_densityProfilePlot.png")

densityProfilePlot(pos, gas_mass, 1UnitfulAstro.Myr, scale = :log10, bins = 50, factor = 6)
savefig(BASE_OUT_PATH * "test_gas_densityProfilePlot.png")

densityProfilePlot(
    [pos, pos],
    [gas_mass, gas_mass],
    1UnitfulAstro.Myr,
    ["sim_1" "sim_2"],
    bins = 50,
    scale = :log10,
    factor = 6,
)
savefig(BASE_OUT_PATH * "test_compare_gas_densityProfilePlot.png")

densityProfilePlot(
    [pos, pos],
    [star_mass, star_mass],
    1UnitfulAstro.Myr,
    ["sim_1" "sim_2"],
    bins = 50,
    scale = :log10,
)
savefig(BASE_OUT_PATH * "test_compare_star_densityProfilePlot.png")

metallicityProfilePlot(pos, gas_mass, gas_z, 1UnitfulAstro.Myr, scale = :log10, bins = 50)
savefig(BASE_OUT_PATH * "test_gas_metallicityProfilePlot.png")

metallicityProfilePlot(pos, star_mass, star_z, 1UnitfulAstro.Myr, scale = :log10, bins = 50)
savefig(BASE_OUT_PATH * "test_star_metallicityProfilePlot.png")

metallicityProfilePlot(
    [pos, pos],
    [gas_mass, gas_mass],
    [gas_z, gas_z],
    1UnitfulAstro.Myr,
    ["sim_1" "sim_2"],
    scale = :log10,
    bins = 50,
)
savefig(BASE_OUT_PATH * "test_compare_gas_metallicityProfilePlot.png")

metallicityProfilePlot(
    [pos, pos],
    [star_mass, star_mass],
    [star_z, star_z],
    1UnitfulAstro.Myr,
    ["sim_1" "sim_2"],
    scale = :log10,
    bins = 50,
)
savefig(BASE_OUT_PATH * "test_compare_star_metallicityProfilePlot.png")

massProfilePlot(pos, star_mass, 1UnitfulAstro.Myr, scale = :log10, bins = 50)
savefig(BASE_OUT_PATH * "test_star_massProfilePlot.png")

massProfilePlot(pos, gas_mass, 1UnitfulAstro.Myr, scale = :log10, bins = 50, factor = 10)
savefig(BASE_OUT_PATH * "test_gas_massProfilePlot.png")

massProfilePlot(
    [pos, pos],
    [gas_mass, gas_mass],
    1UnitfulAstro.Myr,
    ["sim_1" "sim_2"],
    scale = :log10,
    bins = 50,
    factor = 10,
)
savefig(BASE_OUT_PATH * "test_compare_gas_massProfilePlot.png")

massProfilePlot(
    [pos, pos],
    [star_mass, star_mass],
    1UnitfulAstro.Myr,
    ["sim_1" "sim_2"],
    scale = :log10,
    bins = 50,
)
savefig(BASE_OUT_PATH * "test_compare_star_massProfilePlot.png")

sfrTxtPlot(
    sfrtxt_data,
    1,
    [4, 6],
    title = "run_A_01",
    bins = 50,
    scale = (:identity, :log10),
)
savefig(BASE_OUT_PATH * "test_sfrTxtPlot.png")

temperatureHistogramPlot(temp_data, 1UnitfulAstro.Myr, bins = 30)
savefig(BASE_OUT_PATH * "test_temperatureHistogramPlot.png")

rhoTempPlot(temp_data, density, 1UnitfulAstro.Myr)
savefig(BASE_OUT_PATH * "test_rhoTempPlot.png")

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
savefig(BASE_OUT_PATH * "test_KennicuttSchmidtPlot.png")

############################################################################################
# PIPELINE FUNCTIONS.
############################################################################################

gr()

scatterGridPipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "scatter_animation",
    FPS,
    output_path = BASE_OUT_PATH * "scatter_grid/",
    sim_cosmo = SIM_COSMO,
    box_size = BOX_SIZE,
    length_unit = UnitfulAstro.Mpc,
)

densityMapPipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "density_animation",
    FPS,
    output_path = BASE_OUT_PATH * "density_map/",
    sim_cosmo = SIM_COSMO,
    plane = "All",
    box_size = BOX_SIZE,
)

starMapPipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "star_animation",
    FPS,
    output_path = BASE_OUT_PATH * "star_map/",
    sim_cosmo = SIM_COSMO,
    plane = "All",
    box_size = BOX_SIZE,
)

gasStarEvolutionPipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "gas_star_evolution",
    FPS,
    output_path = BASE_OUT_PATH * "gas_star_evolution/",
    sim_cosmo = SIM_COSMO,
    box_size = BOX_SIZE,
)

pgfplotsx()

evolutionSummaryPipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "evolution_summary",
    output_path = BASE_OUT_PATH * "evolution_summary/",
    sim_cosmo = SIM_COSMO,
    mass_factor = 10,
    number_factor = 4,
)

compareSimulationsPipeline(
    [SNAP_NAME, SNAP_NAME],
    [BASE_SRC_PATH, BASE_SRC_PATH],
    ["sim1" "sim2"],
    "compare",
    "star_mass",
    "sfr",
    output_path = BASE_OUT_PATH * "compare_simulations/",
    sim_cosmo = SIM_COSMO,
    title = "SFMS relation",
    x_factor = 10,
)

densityHistogramPipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "density_histogram_animation",
    FPS,
    output_path = BASE_OUT_PATH * "density_histogram/",
    sim_cosmo = SIM_COSMO,
    factor = 10,
)

densityProfilePipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "density_profile_animation",
    FPS,
    "stars";
    output_path = BASE_OUT_PATH * "stars_density_profile/",
    sim_cosmo = SIM_COSMO,
    scale = :log10,
    bins = 80,
    factor = 6,
    box_factor = 0.125,
    box_size = BOX_SIZE,
)

densityProfilePipeline(
    [SNAP_NAME, SNAP_NAME],
    [BASE_SRC_PATH, BASE_SRC_PATH],
    "compare_density_profile_animation",
    FPS,
    "gas",
    ["sim1" "sim2"];
    output_path = BASE_OUT_PATH * "compare_gas_density_profile/",
    sim_cosmo = SIM_COSMO,
    scale = :log10,
    bins = 60,
    factor = 6,
    box_factor = 0.125,
    box_size = BOX_SIZE,
)

metallicityProfilePipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "metallicity_profile_animation",
    FPS,
    "stars";
    output_path = BASE_OUT_PATH * "stars_metallicity_profile/",
    sim_cosmo = SIM_COSMO,
    scale = :log10,
    bins = 80,
    box_factor = 5.0,
    box_size = BOX_SIZE,
)

metallicityProfilePipeline(
    [SNAP_NAME, SNAP_NAME],
    [BASE_SRC_PATH, BASE_SRC_PATH],
    "compare_metallicity_profile_animation",
    FPS,
    "gas",
    ["sim1" "sim2"];
    output_path = BASE_OUT_PATH * "compare_gas_metallicity_profile/",
    sim_cosmo = SIM_COSMO,
    scale = :log10,
    bins = 60,
    box_factor = 5.0,
    box_size = BOX_SIZE,
)

massProfilePipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "mass_profile_animation",
    FPS,
    "stars";
    output_path = BASE_OUT_PATH * "stars_mass_profile/",
    sim_cosmo = SIM_COSMO,
    scale = :log10,
    bins = 80,
    factor = 6,
    box_factor = 0.125,
    box_size = BOX_SIZE,
)

massProfilePipeline(
    [SNAP_NAME, SNAP_NAME],
    [BASE_SRC_PATH, BASE_SRC_PATH],
    "compare_mass_profile_animation",
    FPS,
    "gas",
    ["sim1" "sim2"];
    output_path = BASE_OUT_PATH * "compare_gas_mass_profile/",
    sim_cosmo = SIM_COSMO,
    scale = :log10,
    bins = 60,
    factor = 6,
    box_factor = 0.125,
    box_size = BOX_SIZE,
)

CMDFPipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "CMDF_animation",
    FPS,
    output_path = BASE_OUT_PATH * "CMDF/",
    sim_cosmo = SIM_COSMO,
)

CMDFPipeline(
    [SNAP_NAME, SNAP_NAME],
    [BASE_SRC_PATH, BASE_SRC_PATH],
    "CMDF_animation",
    FPS,
    ["sim1" "sim2"],
    output_path = BASE_OUT_PATH * "compare_CMDF/",
    sim_cosmo = SIM_COSMO,
)

birthHistogramPipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "birth_histogram_animation",
    FPS,
    output_path = BASE_OUT_PATH * "birth_histogram/",
    sim_cosmo = SIM_COSMO,
)

sfrTxtPipeline(
    [SNAP_NAME, SNAP_NAME],
    [BASE_SRC_PATH, BASE_SRC_PATH],
    1,
    [4, 6],
    output_path = BASE_OUT_PATH * "sfr_txt/",
    sim_cosmo = SIM_COSMO,
    title = ["sim_1", "sim_2"],
    bins = 50,
    scale = (:identity, :log10),
)

temperatureHistogramPipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "temperature_histogram_animation",
    FPS,
    output_path = BASE_OUT_PATH * "temperature_histogram/",
    sim_cosmo = SIM_COSMO,
)

rhoTempPipeline(
    SNAP_NAME,
    BASE_SRC_PATH,
    "rho_vs_temp_animation",
    FPS,
    output_path = BASE_OUT_PATH * "rho_vs_temp/",
    sim_cosmo = SIM_COSMO,
)

KennicuttSchmidtPipeline(
    SNAP_NAME,
    BASE_SRC_PATH;
    output_path = BASE_OUT_PATH * "Kennicutt_Schmidt/",
    sim_cosmo = SIM_COSMO,
    box_size = BOX_SIZE,
    bins = 80,
    error_formating = "conf_interval",
    time_unit = UnitfulAstro.yr,
)

############################################################################################
# AUXILIARY FUNCTIONS.
############################################################################################

fig2D = plot(rand(100))
println(GADGETPlotting.relative(fig2D, 0.5, 0.5))
fig3D = surface(rand(100, 100))
println(GADGETPlotting.relative(fig3D, 0.5, 0.5, 0.5))

GADGETPlotting.makeVideo(
    BASE_OUT_PATH * "scatter_grid/images", 
    ".png", 
    BASE_OUT_PATH, 
    "test_video", 
    FPS,
)

x_data = [1:1000...]
y_data = rand(1000)
x_smooth, y_smooth = GADGETPlotting.smoothWindow(x_data, y_data, 50)
plot(x_data, y_data, seriestype = :scatter, legend = false)
plot!(x_smooth, y_smooth, seriestype = :line, lw = 2)
savefig(BASE_OUT_PATH * "test_smoothWindow.png")

positions = pos["gas"]
distances = sqrt.(positions[1, :] .^ 2 .+ positions[2, :] .^ 2 .+ positions[3, :] .^ 2)
box_size = ustrip(Float64, pos["unit"], BOX_SIZE)
r, ρ = GADGETPlotting.densityProfile(gas_mass["mass"], distances, box_size, 80)
plot(r, ρ, lw = 2, xlabel = "r / $(pos["unit"])", ylabel = L"\rho", legend = false)
savefig(BASE_OUT_PATH * "test_densityProfile.png")

r, z = GADGETPlotting.metallicityProfile(
    gas_mass["mass"], 
    distances, 
    gas_z["Z"], 
    box_size, 
    80,
)
plot(r, z, lw = 2, xlabel = "r / $(pos["unit"])", ylabel = "Z / Zsun", legend = false)
savefig(BASE_OUT_PATH * "test_metallicityProfile.png")

r, m = GADGETPlotting.massProfile(gas_mass["mass"], distances, box_size, 80)
plot(r, m, lw = 2, xlabel = "r / $(pos["unit"])", ylabel = "Mass", legend = false)
savefig(BASE_OUT_PATH * "test_massProfile.png")

max_z = findmax(star_z["Z"])
max_Z = max_z[1] / star_mass["mass"][max_z[2]]
z, m = GADGETPlotting.CMDF(star_mass["mass"], star_z["Z"], max_Z, 50)
plot(
    z,
    m,
    lw = 2,
    xlabel = "Z",
    ylabel = L"M_{\star}(< Z) \, / \, M_{\star}",
    legend = false,
)
savefig(BASE_OUT_PATH * "test_CMDF.png")

pos_gas = pos["gas"]
dist_gas = sqrt.(pos_gas[1, :] .^ 2 + pos_gas[2, :] .^ 2)
pos_stars = pos["stars"]
dist_stars = sqrt.(pos_stars[1, :] .^ 2 + pos_stars[2, :] .^ 2)
KSL = GADGETPlotting.KennicuttSchmidtLaw(
    gas_mass["mass"],
    dist_gas,
    temp_data["T"],
    star_mass["mass"],
    dist_stars,
    age_data["ages"],
    ustrip(Float64, temp_data["unit"], 3e4Unitful.K),
    ustrip(Float64, age_data["unit"], 200UnitfulAstro.Myr),
    ustrip(Float64, pos["unit"], BOX_SIZE),
    bins = 80,
)
linear_model = KSL["LM"]
a = round(coef(linear_model)[1], sigdigits = 1)
m = round(coef(linear_model)[2], digits = 1)
a_error = round(stderror(linear_model)[1], sigdigits = 1)
m_error = round(stderror(linear_model)[2], sigdigits = 1)
scatter(KSL["RHO"], KSL["SFR"], label = "Data", xlabel = L"log(\rho)", ylabel = L"log(SFR)")
pl = plot!(KSL["RHO"], predict(linear_model), label = "Fit")
annotate!(
    GADGETPlotting.relative(pl, 0.5, 0.95)..., 
    text(L"SFR = A\,\rho^m", "Courier", 8, :center),
)
annotate!(
    GADGETPlotting.relative(pl, 0.5, 0.9)..., 
    text(L"m = %$m \pm %$m_error", "Courier", 8, :center),
)
annotate!(
    GADGETPlotting.relative(pl, 0.5, 0.85)...,
    text(L"log(A) = %$a \pm %$a_error", "Courier", 8, :center),
)
savefig(BASE_OUT_PATH * "test_KennicuttSchmidtLaw.png")

println(GADGETPlotting.format_error(69.42069, 0.038796))
println(GADGETPlotting.format_error(69.42069, 0.018796))
println(GADGETPlotting.format_error(69.42069, 0.0))
println(GADGETPlotting.format_error(69.42069, 73.4))

############################################################################################
# DELETE ALL TESTING FILES PRODUCED.
############################################################################################

# rm(BASE_OUT_PATH, recursive = true)

println("Everything worked just fine!!")