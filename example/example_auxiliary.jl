############################################################################################
# Example auxiliary functions
############################################################################################

##############
## Data
##############

pos_gas = pos["gas"]
pos_stars = pos["stars"]
m_gas = gas_mass["mass"]
dist_gas_2D_3D = sqrt.(pos_gas[1, :] .^ 2 .+ pos_gas[2, :] .^ 2 .+ pos_gas[3, :] .^ 2)
dist_gas_2D = sqrt.(pos_gas[1, :] .^ 2 + pos_gas[2, :] .^ 2)
dist_stars = sqrt.(pos_stars[1, :] .^ 2 + pos_stars[2, :] .^ 2)
box_size = ustrip(Float64, pos["unit"], BOX_SIZE)
max_z = findmax(star_z["Z"])
max_Z = max_z[1] / star_mass["mass"][max_z[2]]

##############
## Functions
##############

relative_2D = GADGETPlotting.relative(plot(rand(100)), 0.5, 0.5)
relative_3D = GADGETPlotting.relative(surface(rand(100, 100)), 0.5, 0.5, 0.5)

# GADGETPlotting.make_video(
#     joinpath(BASE_OUT_PATH, "scatter_grid/images"),
#     ".png",
#     BASE_OUT_PATH,
#     "test_video",
#     FPS,
# )

smooth_w = GADGETPlotting.smooth_window([1:1000...], rand(1000), 50)

box_size = ustrip(Float64, pos["unit"], BOX_SIZE)

density_p = GADGETPlotting.density_profile(m_gas, dist_gas_2D_3D, box_size, 80)

metallicity_p = GADGETPlotting.metallicity_profile(
    m_gas, 
    dist_gas_2D_3D, 
    gas_z["Z"], 
    box_size, 
    80,
)

mass_p = GADGETPlotting.mass_profile(m_gas, dist_gas_2D_3D, box_size, 80)

cmdf = GADGETPlotting.compute_cmdf(star_mass["mass"], star_z["Z"], max_Z, 50)

ksl = GADGETPlotting.kennicutt_schmidt_law(
    m_gas,
    dist_gas_2D,
    temp_data["temperature"],
    star_mass["mass"],
    dist_stars,
    age_data["ages"],
    ustrip(Float64, temp_data["unit"], 3e4Unitful.K),
    ustrip(Float64, age_data["unit"], 200UnitfulAstro.Myr),
    ustrip(Float64, pos["unit"], BOX_SIZE),
    bins = 80,
)

format_err_1 = GADGETPlotting.format_error(69.42069, 0.038796)
format_err_2 = GADGETPlotting.format_error(69.42069, 0.018796)
format_err_3 = GADGETPlotting.format_error(69.42069, 0.0)
format_err_4 = GADGETPlotting.format_error(69.42069, 73.4)

pass = GADGETPlotting.pass_all(FIRST_SNAP, "gas")

energy_i = GADGETPlotting.energy_integrand(read_header(FIRST_SNAP), 1.0)

num_int_1 = GADGETPlotting.num_integrate(sin, 0, 3π)
num_int_2 = GADGETPlotting.num_integrate(x -> x^3 + 6 * x^2 + 9 * x + 2, 0, 4.69)
num_int_3 = GADGETPlotting.num_integrate(x -> exp(x^x), 0, 1)
num_int_4 = GADGETPlotting.num_integrate(x -> sqrt(sqrt(1 / (x + 1))), 0, 1)

rc = GADGETPlotting.center_of_mass(pos_gas, m_gas)

##############
## Testing 
##############

jldsave(
    joinpath(BASE_OUT_PATH, "data_auxiliary.jld2"); 
    relative_2D, relative_3D, smooth_w,
    density_p, metallicity_p, mass_p,
    cmdf, ksl, format_err_1,
    format_err_2, format_err_3, format_err_4,
    pass, energy_i, num_int_1, 
    num_int_2, num_int_3, num_int_4, rc,
)