########################################################################################
# Testing script for GADGETPlotting.jl, it can be run as is, and it shouldn't 
# throw any errors.
########################################################################################

include("GADGETPlotting.jl")

"Base path for the directories where the figures and animations will be saved."
const BASE_OUT_PATH = "results/"  

"Directory containing the snapshot files."   
const SIM_PATH = "../../sim_data/"  

const SNAP_PATH = ["run_A_01/", "run_C_01/", "run_D_01/", "run_E_01/"]

"Base name of the snapshot files."
const SNAP_NAME = ["snap", "snap", "snap", "snap"] 

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

snap_files = []
for simu in SNAP_PATH
    snaps = getSnapshots(SNAP_NAME[1], SIM_PATH * simu)
    push!(snap_files, snaps["snap_files"])
end

pos = Dict{String, Any}[]
gas_mass = Dict{String, Any}[]
gas_z = Dict{String, Any}[]
for file in snap_files
    push!(pos, positionData(file[151], sim_cosmo=SIM_COSMO, length_unit=UnitfulAstro.kpc, box_size=BOX_SIZE))
    push!(gas_mass, massData(file[151], "gas", sim_cosmo=SIM_COSMO))
    push!(gas_z, zData(file[151], "gas", sim_cosmo=SIM_COSMO))
end

densityProfilePlot(pos, gas_mass, 1 * UnitfulAstro.Myr, ["run_A_01" "run_C_01" "run_D_01" "run_E_01"], scale=:log10, bins=60, factor=8, region_factor=0.1)
Base.invokelatest(savefig, BASE_OUT_PATH * "comparison.png")

metallicityProfilePlot(pos, gas_mass, gas_z, 1 * UnitfulAstro.Myr, ["run_A_01" "run_C_01" "run_D_01" "run_E_01"], scale=:log10, bins=60, region_factor=0.3)
Base.invokelatest(savefig, BASE_OUT_PATH * "comparison_2.png")