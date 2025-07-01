push!(LOAD_PATH, "./src/")
using Documenter, GADGETPlotting
using DocumenterTools: Themes

# True if being deployed, false if being compile locally
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

# Compile documentation
makedocs(
    sitename="GADGETPlotting.jl",
    authors = "Ezequiel Lozano",
    format = Documenter.HTML(
        assets=["./assets/GI-docs.css"],
        prettyurls = CI,
        edit_link = "main",
    ),
    modules = [GADGETPlotting],
    pages = [
        "Introduction"                  => "index.md",
        "Pipeline Functions"            => "pipelines.md",
        "Plotting Functions"            => "plotting.md",
        "Data Acquisition Functions"    => "data_acquisition.md",
        "Auxiliary Functions"           => "auxiliary.md",
        "Index"                         => "func_list.md",
    ],
)

# Deploy documentation to GitHub pages
if CI
    deploydocs(
        repo = "github.com/Ezequiel92/GADGETPlotting.git",
        devbranch = "main",
        versions = ["dev" => "dev"],
    )
end
