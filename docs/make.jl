push!(LOAD_PATH, "./src/")
using Documenter, GADGETPlotting

const formats = Any[
    Documenter.HTML(
        edit_branch = "main",
    ),
]

makedocs(
    sitename="GADGETPlotting.jl",
    format = Documenter.HTML(),
    modules = [GADGETPlotting],
)

deploydocs(
    repo = "github.com/Ezequiel92/GADGETPlotting.git",
    devbranch = "main",
	devurl = "stable",
)
