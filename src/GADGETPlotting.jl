############################################################################################
# Julia module for creating plots, GIFs, and videos from the data 
# produced by GAGET2/3/4 simulations. 
############################################################################################

module GADGETPlotting

using GadgetIO, GadgetUnits, SPHtoGrid, SPHKernels
using Unitful, UnitfulAstro
using Plots, LaTeXStrings, StatsPlots.PlotMeasures, AverageShiftedHistograms, GLM
using Glob, FileIO, VideoIO, DelimitedFiles, Accessors, ProgressMeter, LinearAlgebra

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
    @eval Base.Experimental.@optlevel 3
end

############################################################################################
# Global constants. 
############################################################################################

@doc raw"""
``H_0 = 100 \, \mathrm{km} \, \mathrm{s}^{-1} \, \mathrm{Mpc}^{-1} \ \mathrm{in} \ \mathrm{Gyr}^{-1}``
"""
const HUBBLE_CONST = 0.102201

"""
Solar metallicity.

M. Asplund et al. (2009). *The Chemical Composition of the Sun.* Annual Review of Astronomy 
and Astrophysics, **47(1)**, 481–522. [https://doi.org/10.1146/annurev.astro.46.060407.145222](https://doi.org/10.1146/annurev.astro.46.060407.145222)
"""
const SOLAR_METALLICITY = 0.0134

"""
Slope for the Kennicutt-Schmidt law.

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies.* The Astrophysical 
Journal, **498(2)**, 541-552. [https://doi.org/10.1086/305588](https://doi.org/10.1086/305588)
"""
const KENNICUTT98_SLOPE = 1.4

"""
Intercept for the Kennicutt-Schmidt law.

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies.* The Astrophysical 
Journal, **498(2)**, 541-552. [https://doi.org/10.1086/305588](https://doi.org/10.1086/305588)
"""
const KENNICUTT98_INTERCEPT = 2.5e-4 * (UnitfulAstro.Msun / UnitfulAstro.yr / UnitfulAstro.kpc^2)

"""
Unit of area density for the Kennicutt-Schmidt law.

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies.* The Astrophysical 
Journal, **498(2)**, 541-552. [https://doi.org/10.1086/305588](https://doi.org/10.1086/305588)
"""
const KENNICUTT98_RHO_UNIT = 1.0 * UnitfulAstro.Msun / UnitfulAstro.pc^2

############################################################################################
# Functions. 
############################################################################################

include("auxiliary.jl")
include("data_acquisition.jl")
include("plotting.jl")
include("pipelines.jl")
 
export 
    # Plotting functions ###################################################################
    scatter_grid_plot,                              
    density_map_plot,
    star_map_plot,
    gas_star_evolution_plot,
    cmdf_plot,
    birth_histogram_plot,
    time_series_plot,
    scale_factor_series_plot,
    redshift_series_plot,
    compare_simulations_plot,
    density_histogram_plot,
    density_profile_plot,
    metallicity_profile_plot,
    mass_profile_plot,
    sfr_txt_plot,
    temperature_histogram_plot,
    rho_temp_plot,
    kennicutt_schmidt_plot,  
    cpu_txt_plot,  
    quantities_2D_plot,           
    # Pipeline functions ###################################################################           
    scatter_grid_pipeline,                
    density_map_pipeline,
    star_map_pipeline,
    gas_star_evolution_pipeline,
    evolution_summary_pipeline,
    compare_simulations_pipeline,
    density_histogram_pipeline,
    density_profile_pipeline,
    metallicity_profile_pipeline,
    mass_profile_pipeline,
    cmdf_pipeline,
    birth_histogram_pipeline,
    sfr_txt_pipeline,
    temperature_histogram_pipeline,
    rho_temp_pipeline,
    kennicutt_schmidt_pipeline,   
    cpu_txt_pipeline,   
    quantities_2D_pipeline,     
    # Auxiliary functions ##################################################################
    comparison,                     
    deep_comparison

end