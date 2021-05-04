############################################################################################
# AUXILIARY FUNCTIONS
############################################################################################

relative_2D = GADGETPlotting.relative(plot(rand(100)), 0.5, 0.5)
relative_3D = GADGETPlotting.relative(surface(rand(100, 100)), 0.5, 0.5, 0.5)

# GADGETPlotting.makeVideo(
#     joinpath(BASE_OUT_PATH, "scatter_grid/images"),
#     ".png",
#     BASE_OUT_PATH,
#     "test_video",
#     FPS,
# )

# pgfplotsx()

smooth_w = GADGETPlotting.smoothWindow([1:1000...], rand(1000), 50)

positions = pos["gas"]
distances = sqrt.(positions[1, :] .^ 2 .+ positions[2, :] .^ 2 .+ positions[3, :] .^ 2)
box_size = ustrip(Float64, pos["unit"], BOX_SIZE)
density_p = GADGETPlotting.densityProfile(gas_mass["mass"], distances, box_size, 80)

metallicity_p = GADGETPlotting.metallicityProfile(
    gas_mass["mass"], 
    distances, 
    gas_z["Z"], 
    box_size, 
    80,
)

mass_p = GADGETPlotting.massProfile(gas_mass["mass"], distances, box_size, 80)

max_z = findmax(star_z["Z"])
max_Z = max_z[1] / star_mass["mass"][max_z[2]]
cmdf = GADGETPlotting.CMDF(star_mass["mass"], star_z["Z"], max_Z, 50)

pos_gas = pos["gas"]
dist_gas = sqrt.(pos_gas[1, :] .^ 2 + pos_gas[2, :] .^ 2)
pos_stars = pos["stars"]
dist_stars = sqrt.(pos_stars[1, :] .^ 2 + pos_stars[2, :] .^ 2)
ksl = GADGETPlotting.KennicuttSchmidtLaw(
    gas_mass["mass"],
    dist_gas,
    temp_data["temperature"],
    star_mass["mass"],
    dist_stars,
    age_data["ages"],
    ustrip(Float64, temp_data["unit"], 3e4Unitful.K),
    ustrip(Float64, age_data["unit"], 200UnitfulAstro.Myr),
    ustrip(Float64, pos["unit"], BOX_SIZE),
    bins = 80,
)
# linear_model = KSL["LM"]
# a = round(coef(linear_model)[1], sigdigits = 1)
# m = round(coef(linear_model)[2], digits = 1)
# a_error = round(stderror(linear_model)[1], sigdigits = 1)
# m_error = round(stderror(linear_model)[2], sigdigits = 1)
# scatter(KSL["RHO"], KSL["SFR"], label = "Data", xlabel = L"log(\rho)", ylabel = L"log(SFR)")
# pl = plot!(KSL["RHO"], predict(linear_model), label = "Fit")
# annotate!(
#     GADGETPlotting.relative(pl, 0.5, 0.95)...,
#     text(L"SFR = A\,\rho^m", "Courier", 8, :center),
# )
# annotate!(
#     GADGETPlotting.relative(pl, 0.5, 0.9)...,
#     text(L"m = %$m \pm %$m_error", "Courier", 8, :center),
# )
# annotate!(
#     GADGETPlotting.relative(pl, 0.5, 0.85)...,
#     text(L"log(A) = %$a \pm %$a_error", "Courier", 8, :center),
# )
# savefig(joinpath(BASE_OUT_PATH, "test_KennicuttSchmidtLaw.png"))

format_err_1 = GADGETPlotting.format_error(69.42069, 0.038796)
format_err_2 = GADGETPlotting.format_error(69.42069, 0.018796)
format_err_3 = GADGETPlotting.format_error(69.42069, 0.0)
format_err_4 = GADGETPlotting.format_error(69.42069, 73.4)

pass = GADGETPlotting.pass_all(joinpath(BASE_SRC_PATH, SNAP_NAME * "_000"), "gas")

header = read_header(joinpath(BASE_SRC_PATH, SNAP_NAME * "_000"))
energy_i = GADGETPlotting.energy_integrand(header, 1.0)

num_int_1 = GADGETPlotting.num_integrate(sin, 0, 3π)
num_int_2 = GADGETPlotting.num_integrate(x -> x^3 + 6 * x^2 + 9 * x + 2, 0, 4.69)
num_int_3 = GADGETPlotting.num_integrate(x -> exp(x^x), 0, 1)
num_int_4 = GADGETPlotting.num_integrate(x -> sqrt(sqrt(1 / (x + 1))), 0, 1)

# Save numerical result for tests
jldsave(
    joinpath(BASE_OUT_PATH, "data_auxiliary.jld2"); 
    relative_2D, relative_3D, smooth_w,
    density_p, metallicity_p, mass_p,
    cmdf, ksl, format_err_1,
    format_err_2, format_err_3, format_err_4,
    pass, energy_i, num_int_1, 
    num_int_2, num_int_3, num_int_4,
)