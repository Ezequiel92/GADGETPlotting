using Plots: display
############################################################################################
#   Example scripts for GADGETPlotting.jl
#
# - When run as is, it shouldn't throw any errors
# - This script shows how to import GADGETPlotting.jl and gives examples on how 
#   to use every function, for an isolated and a cosmological simulation
# - When run it produces `.jld2` files, used for testing.
############################################################################################

push!(LOAD_PATH, joinpath(@__DIR__, "../src/"))
using GADGETPlotting, GadgetIO, Unitful, UnitfulAstro, Plots, LaTeXStrings, GLM, JLD2

"Base name of the snapshot files, set in the GADGET variable SnapshotFileBase."
const SNAP_NAME = "snap"

############################################################################################
# Examples (isolated)
############################################################################################

"Base path to the directories where the output images and animations will be saved."
BASE_OUT_PATH = joinpath(@__DIR__, "example_results/isolated")

"Directory containing the snapshot files."
BASE_SRC_PATH = joinpath(@__DIR__, "example_data/isolated")

"Path to the fist snapshot."
FIRST_SNAP = joinpath(BASE_SRC_PATH, SNAP_NAME * "_000")

"Side dimension of the simulated region, for the case of vacuum boundary conditions."
BOX_SIZE = 200.0UnitfulAstro.kpc

"Value of ComovingIntegrationOn: 0 -> Newtonian simulation, 1 -> Cosmological simulation."
SIM_COSMO = 0

"Frame rate for the animations."
FPS = 4

"One particular snapshot index for testing purposes."
SNAP_N = 21

mkpath(BASE_OUT_PATH)

include("example_data_acquisition.jl")
include("example_plotting.jl")
include("example_pipeline.jl")
include("example_auxiliary.jl")

############################################################################################
# Examples (cosmological)
############################################################################################

BASE_OUT_PATH = joinpath(@__DIR__, "example_results/cosmological")
BASE_SRC_PATH = joinpath(@__DIR__, "example_data/cosmological")
FIRST_SNAP = joinpath(BASE_SRC_PATH, SNAP_NAME * "dir_019/" * SNAP_NAME * "_019")
SIM_COSMO = 1
FPS = 1
SNAP_N = 2 

mkpath(BASE_OUT_PATH)

include("example_data_acquisition.jl")
include("example_plotting.jl")
include("example_pipeline.jl")
include("example_auxiliary.jl")

println("Everything worked just fine!!")
