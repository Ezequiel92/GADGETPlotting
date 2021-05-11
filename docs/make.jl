push!(LOAD_PATH, "./src/")
using Documenter, GADGETPlotting

makedocs(
    sitename="GADGETPlotting.jl",
    format = Documenter.HTML(),
    modules = [GADGETPlotting],
)

deploydocs(
    repo = "github.com/Ezequiel92/GADGETPlotting.jl.git",
	devbranch = "main",
)