push!(LOAD_PATH, "./src/")
using GADGETPlotting, Unitful, UnitfulAstro, Plots, LaTeXStrings, GLM

"Base path to the directories where the output images and animations will be saved."
const BASE_OUT_PATH = joinpath(@__DIR__, "example_results")

"Directory containing the snapshot files."
const BASE_SRC_PATH = joinpath(@__DIR__, "test_data")

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
# Examples.
############################################################################################

include("example_data_acquisition.jl")
include("example_plotting.jl")
include("example_pipeline.jl")
include("example_auxiliary.jl")

println("Everything worked just fine!!")

############################################################################################
# Delete all testing files produced.
############################################################################################

# rm(BASE_OUT_PATH, recursive = true)
