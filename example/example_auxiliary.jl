############################################################################################
# AUXILIARY FUNCTIONS.
############################################################################################

fig2D = plot(rand(100))
println(GADGETPlotting.relative(fig2D, 0.5, 0.5))
fig3D = surface(rand(100, 100))
println(GADGETPlotting.relative(fig3D, 0.5, 0.5, 0.5))

GADGETPlotting.makeVideo(
    joinpath(BASE_OUT_PATH, "scatter_grid/images"),
    ".png",
    BASE_OUT_PATH,
    "test_video",
    FPS,
)

pgfplotsx()

x_data = [1:1000...]
y_data = rand(1000)
x_smooth, y_smooth = GADGETPlotting.smoothWindow(x_data, y_data, 50)
plot(x_data, y_data, seriestype = :scatter, legend = false)
plot!(x_smooth, y_smooth, seriestype = :line, lw = 2)
savefig(joinpath(BASE_OUT_PATH, "test_smoothWindow.png"))

positions = pos["gas"]
distances = sqrt.(positions[1, :] .^ 2 .+ positions[2, :] .^ 2 .+ positions[3, :] .^ 2)
box_size = ustrip(Float64, pos["unit"], BOX_SIZE)
r, ρ = GADGETPlotting.densityProfile(gas_mass["mass"], distances, box_size, 80)
plot(r, ρ, lw = 2, xlabel = "r / $(pos["unit"])", ylabel = L"\rho", legend = false)
savefig(joinpath(BASE_OUT_PATH, "test_densityProfile.png"))

r, z = GADGETPlotting.metallicityProfile(
    gas_mass["mass"], 
    distances, 
    gas_z["Z"], 
    box_size, 
    80,
)
plot(r, z, lw = 2, xlabel = "r / $(pos["unit"])", ylabel = "Z / Zsun", legend = false)
savefig(joinpath(BASE_OUT_PATH, "test_metallicityProfile.png"))

r, m = GADGETPlotting.massProfile(gas_mass["mass"], distances, box_size, 80)
plot(r, m, lw = 2, xlabel = "r / $(pos["unit"])", ylabel = "Mass", legend = false)
savefig(joinpath(BASE_OUT_PATH, "test_massProfile.png"))

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
savefig(joinpath(BASE_OUT_PATH, "test_CMDF.png"))

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
savefig(joinpath(BASE_OUT_PATH, "test_KennicuttSchmidtLaw.png"))

println(GADGETPlotting.format_error(69.42069, 0.038796))
println(GADGETPlotting.format_error(69.42069, 0.018796))
println(GADGETPlotting.format_error(69.42069, 0.0))
println(GADGETPlotting.format_error(69.42069, 73.4))

display(GADGETPlotting.pass_all(joinpath(BASE_SRC_PATH, SNAP_NAME * "_000"), "gas"))

header = read_header(joinpath(BASE_SRC_PATH, SNAP_NAME * "_000"))
println(GADGETPlotting.energy_integrand(header, 0.1))
println(GADGETPlotting.energy_integrand(header, 0.5))
println(GADGETPlotting.energy_integrand(header, 0.8))

println(GADGETPlotting.num_integrate(sin, 0, 3π))
println(GADGETPlotting.num_integrate(x -> x^3 + 6 * x^2 + 9 * x + 2, 0, 4.69))
println(GADGETPlotting.num_integrate(x -> exp(x^x), 0, 1))