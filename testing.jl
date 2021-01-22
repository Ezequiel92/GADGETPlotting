########################################################################################
# Testing script for GADGETPlotting.jl, it can be run as is, and it shouldn't 
# throw any errors.
########################################################################################

include("GADGETPlotting.jl")

"Base path for the directories where the figures and animations will be saved."
const BASE_OUT_PATH = "results/"  

"Directory containing the snapshot files."   
const SNAP_PATH = "test_snapshots/"  

"Base name of the snapshot files."
const SNAP_NAME = "snap" 

"Side dimension of the simulated region, with units, for the case of vacuum boundary conditions."
const BOX_SIZE = 200 * UnitfulAstro.kpc 

"Value of ComovingIntegrationOn: 0 -> Newtonian simulation, 1 -> Cosmological simulation."               
const SIM_COSMO = 0 

"Frame rate for the animations."                   
const FPS = 4    

mkpath(BASE_OUT_PATH)

########################################################################################
# TEST OF DATA ACQUISITION FUNCTIONS.
########################################################################################

snaps = getSnapshots(SNAP_NAME, SNAP_PATH)
snap_files = snaps["snap_files"]
display(snaps)
println()

time_series = timeSeriesData(snap_files, sim_cosmo=SIM_COSMO)
display(time_series)
println()

pos = positionData( snap_files[21], 
                    sim_cosmo=SIM_COSMO, 
                    length_unit=UnitfulAstro.Mpc, 
                    box_size=BOX_SIZE)
display(pos)
println()

density = densityData(snap_files[21], sim_cosmo=SIM_COSMO)
display(density)
println()

hsml = hsmlData(snap_files[21], sim_cosmo=SIM_COSMO)
display(hsml)
println()

gas_mass = massData(snap_files[21], "gas", sim_cosmo=SIM_COSMO)
display(gas_mass)
println()

dm_mass = massData(snap_files[21], "dark_matter", sim_cosmo=SIM_COSMO)
display(dm_mass)
println()

star_mass = massData(snap_files[21], "stars", sim_cosmo=SIM_COSMO)
display(star_mass)
println()

gas_z = zData(snap_files[21], "gas", sim_cosmo=SIM_COSMO)
display(gas_z)
println()

sfrtxtdata = sfrTxtData(SNAP_PATH, SNAP_NAME, sim_cosmo=SIM_COSMO)
display(sfrtxtdata)

########################################################################################
# TEST OF PLOTTING FUNCTIONS.
########################################################################################

scatterGridPlot(pos)
savefig(BASE_OUT_PATH * "test_scatterGridPlot.png")

densityMapPlot(pos, gas_mass, hsml, density, plane="XY")
savefig(BASE_OUT_PATH * "test_densityMapPlot_XY.png")

densityMapPlot(pos, gas_mass, hsml, density, plane="All")
savefig(BASE_OUT_PATH * "test_densityMapPlot_All.png")

color = [:batlow, :bone, :CMRmap, :grayC, :inferno, :seaborn_rocket_gradient, :YlOrRd_9]

for c in color
    densityMapPlot(pos, gas_mass, hsml, density, plane="All", color=c)
    savefig(BASE_OUT_PATH * "test_densityMapPlot_" * string(c) * ".png")
end

starMapPlot(pos, plane="XY")
savefig(BASE_OUT_PATH * "test_starMapPlot_XY.png")

starMapPlot(pos, plane="XZ")
savefig(BASE_OUT_PATH * "test_starMapPlot_XZ.png")

starMapPlot(pos)
savefig(BASE_OUT_PATH * "test_starMapPlot_All.png")

gasStarEvolutionPlot(time_series, pos, 21)
savefig(BASE_OUT_PATH * "test_gasStarEvolutionPlot.png")

timeSeriesPlot(time_series, mass_factor=0, number_factor=4)
Base.invokelatest(savefig, BASE_OUT_PATH * "test_timeSeriesPlot.png")

scaleFactorSeriesPlot(time_series, mass_factor=0, number_factor=4)
Base.invokelatest(savefig, BASE_OUT_PATH * "test_scaleFactorSeriesPlot.png")

redshiftSeriesPlot(time_series, mass_factor=0, number_factor=4)
Base.invokelatest(savefig, BASE_OUT_PATH * "test_redshiftSeriesPlot.png")

compareSimulationsPlot( "star_mass", 
                        "sfr", 
                        ["sim1" "sim2"], 
                        (time_series, time_series); 
                        title="SFMS relation", 
                        x_factor=10,
                        log_scale=true) 
Base.invokelatest(savefig, BASE_OUT_PATH * "test_compareSimulationsPlot.png")

densityHistogramPlot(density, 1 * UnitfulAstro.Myr, factor=10)
Base.invokelatest(savefig, BASE_OUT_PATH * "test_densityHistogramPlot.png")

densityProfilePlot(pos, gas_mass,  1 * UnitfulAstro.Myr, bins=50, factor=13)
Base.invokelatest(savefig, BASE_OUT_PATH * "test_densityProfilePlot.png")

metallicityProfilePlot(pos, gas_mass, gas_z, 1 * UnitfulAstro.Myr, bins=100)
Base.invokelatest(savefig, BASE_OUT_PATH * "test_metallicityProfilePlot.png")

sfrTxtPlot( SNAP_PATH, 
            SNAP_NAME, 
            1, 
            [4, 6], 
            title="run_A_01", 
            bins=50, 
            scale=[:identity, :log10], 
            sim_cosmo=SIM_COSMO)
Base.invokelatest(savefig, BASE_OUT_PATH * "test_sfrTxtPlot.png")

########################################################################################
# TEST OF PIPELINE FUNCTIONS.
########################################################################################

scatterGridPipeline(SNAP_NAME, 
                    SNAP_PATH, 
                    BASE_OUT_PATH, 
                    "scatter_anim", 
                    FPS,
                    sim_cosmo=SIM_COSMO,
                    length_unit=UnitfulAstro.Mpc,
                    region_size=BOX_SIZE)

densityMapPipeline( SNAP_NAME, 
                    SNAP_PATH, 
                    BASE_OUT_PATH, 
                    "density_anim", 
                    FPS,
                    plane="All",
                    sim_cosmo=SIM_COSMO,
                    region_size=BOX_SIZE)

starMapPipeline(SNAP_NAME, 
                SNAP_PATH, 
                BASE_OUT_PATH, 
                "star_anim", 
                FPS,
                plane="All",
                sim_cosmo=SIM_COSMO,
                region_size=BOX_SIZE)

gasStarEvolutionPipeline(   SNAP_NAME,
                            SNAP_PATH, 
                            BASE_OUT_PATH,
                            "gas_star_evolution",
                            FPS,
                            sim_cosmo=SIM_COSMO,
                            region_size=BOX_SIZE)

pgfplotsx()

evolutionSummaryPipeline(   SNAP_NAME, 
                            SNAP_PATH,  
                            BASE_OUT_PATH, 
                            "evolution_summary",
                            sim_cosmo=SIM_COSMO,
                            mass_factor=10, 
                            number_factor=4)

compareSimulationsPipeline( [SNAP_NAME, SNAP_NAME],
                            [SNAP_PATH, SNAP_PATH],
                            BASE_OUT_PATH, 
                            ["sim1" "sim2"], 
                            "compare",
                            "star_mass",
                            "sfr",
                            x_factor=10,
                            folder="compare_simulations/",
                            title="SFMS relation")

densityHistogramPipeline(   SNAP_NAME, 
                            SNAP_PATH, 
                            BASE_OUT_PATH, 
                            "density_histogram_anim", 
                            FPS,
                            sim_cosmo=SIM_COSMO,
                            factor=10)

densityProfilePipeline( SNAP_NAME, 
                        SNAP_PATH, 
                        BASE_OUT_PATH, 
                        "density_profile_anim", 
                        FPS,
                        "gas";
                        bins=80,
                        factor=6,
                        region_factor=0.125,
                        sim_cosmo=SIM_COSMO,
                        region_size=BOX_SIZE)

metallicityProfilePipeline( SNAP_NAME, 
                            SNAP_PATH, 
                            BASE_OUT_PATH, 
                            "metallicity_profile_anim", 
                            FPS,
                            "gas";
                            bins=80,
                            region_factor=5.0,
                            sim_cosmo=SIM_COSMO,
                            region_size=BOX_SIZE)

########################################################################################
# TEST OF AUXILIARY FUNCTIONS.
########################################################################################

println(uglyUnitDeserialization(UnitfulAstro.Msun))
println(uglyUnitDeserialization(UnitfulAstro.Myr))
println(uglyUnitDeserialization(UnitfulAstro.kpc))

fig2D = plot(rand(100))
println(relative(fig2D, 0.5, 0.5))
fig3D = surface(rand(100, 100))
println(relative(fig3D, 0.5, 0.5, 0.5))

makeVideo(BASE_OUT_PATH * "scatter_grid/images", ".png", BASE_OUT_PATH, "test_video", FPS)

x_data = [1:1000...]
y_data = rand(1000)
x_smooth, y_smooth = smoothWindow(x_data, y_data, 50)
plot(x_data, y_data, seriestype=:scatter, legend=false)
plot!(x_smooth, y_smooth, seriestype=:line, lw=2)
savefig(BASE_OUT_PATH * "test_smoothWindow.png")

positions = pos["gas"]
distances = sqrt.(positions[:,1].^2 .+ positions[:,2].^2 .+  positions[:,3].^2)
box_size = ustrip(Float64, unit(BOX_SIZE), BOX_SIZE)
r, ρ = densityProfile(gas_mass["mass"], distances, box_size, 80)
plot(r, ρ, lw=2, xlabel="r / $(pos["unit"])", ylabel=L"\rho", legend=false)
savefig(BASE_OUT_PATH * "test_densityProfile.png")

positions = pos["gas"]
distances = sqrt.(positions[:,1].^2 .+ positions[:,2].^2 .+  positions[:,3].^2)
box_size = ustrip(Float64, unit(BOX_SIZE), BOX_SIZE)
r, z = metallicityProfile(gas_mass["mass"], distances, gas_z["Z"], box_size, 80)
plot(r, z, lw=2, xlabel="r / $(pos["unit"])", ylabel="Z / Zsun", legend=false)
savefig(BASE_OUT_PATH * "test_metallicityProfile.png")


println("Everything worked just fine!!")


########################################################################################
# DELETE ALL GENERATED TESTING FILES.
########################################################################################

# rm(BASE_OUT_PATH, recursive=true)

########################################################################################
# EXTRA TESTS.
# This functions when executed should throw an error.
########################################################################################

# massData(snap_files[21], "unobtainium", sim_cosmo=SIM_COSMO)
# zData(snap_files[21], "unobtainium", sim_cosmo=SIM_COSMO)