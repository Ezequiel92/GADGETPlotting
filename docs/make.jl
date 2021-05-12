push!(LOAD_PATH, "./src/")
using Documenter, GADGETPlotting
using DocumenterTools: Themes

CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

# Themes.compile(joinpath(@__DIR__,"src/assets/geostats-light.scss"), joinpath(@__DIR__,"src/assets/themes/documenter-light.css"))
# Themes.compile(joinpath(@__DIR__,"src/assets/geostats-dark.scss"), joinpath(@__DIR__,"src/assets/themes/documenter-dark.css"))

for file in ("geostats-light.scss", "geostats-dark.scss")
    download("https://raw.githubusercontent.com/JuliaEarth/GeoStats.jl/master/docs/src/assets/$file", joinpath(@__DIR__, "src/assets/$file"))
end

# compile the themes
Themes.compile(joinpath(@__DIR__, "src/assets/geostats-light.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-light.css"))
Themes.compile(joinpath(@__DIR__, "src/assets/geostats-dark.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"))

#isdir(datadir()) && rm(datadir(); force = true, recursive = true)

makedocs(
    sitename="GADGETPlotting.jl",
    authors = "Ezequiel Lozano",
    format = Documenter.HTML(
        assets=[asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css)],
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

if CI
    deploydocs(
        repo = "github.com/Ezequiel92/GADGETPlotting.git",
        devbranch = "main",
        versions = ["dev" => "dev"],
    ) 
end


