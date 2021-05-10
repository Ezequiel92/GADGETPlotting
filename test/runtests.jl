push!(LOAD_PATH, joinpath(@__DIR__, "../src/"))
using GADGETPlotting, GadgetIO, Test, ReferenceTests, Plots, JLD2, UnitfulAstro, Unitful

const SNAP_NAME = "snap"

BASE_SRC_PATH = joinpath(@__DIR__, "../example/example_data/isolated")
BASE_DATA_PATH = joinpath(@__DIR__, "../example/example_results/isolated")
FIRST_SNAP = joinpath(BASE_SRC_PATH, "snap_000")
BOX_SIZE = 200UnitfulAstro.kpc
SIM_COSMO = 0
FPS = 4
SNAP_N = 21

println("Testing with an isolated simulation...\n")

include("test_data_acquisition.jl")
include("test_plotting.jl")
include("test_auxiliary.jl")

BASE_SRC_PATH = joinpath(@__DIR__, "../example/example_data/cosmological")
BASE_DATA_PATH = joinpath(@__DIR__, "../example/example_results/cosmological")
FIRST_SNAP = joinpath(BASE_SRC_PATH, SNAP_NAME * "dir_019/" * SNAP_NAME * "_019")
SIM_COSMO = 1
FPS = 1
SNAP_N = 2 

println("\nTesting with a cosmological simulation...\n")

include("test_data_acquisition.jl")
include("test_plotting.jl")
include("test_auxiliary.jl")