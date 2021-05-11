# GADGETPlotting.jl

```@contents
Depth = 3
```
## Exported Functions

### Pipelines

```@docs
scatter_grid_pipeline               
density_map_pipeline
star_map_pipeline
gas_star_evolution_pipeline
evolution_summary_pipeline
compare_simulations_pipeline
density_histogram_pipeline
density_profile_pipeline
metallicity_profile_pipeline
mass_profile_pipeline
cmdf_pipeline
birth_histogram_pipeline
sfr_txt_pipeline
temperature_histogram_pipeline
rho_temp_pipeline
kennicutt_schmidt_pipeline
```

### Plotting

```@docs
scatter_grid_plot                             
density_map_plot
star_map_plot
gas_star_evolution_plot
cmdf_plot
birth_histogram_plot
time_series_plot
scale_factor_series_plot
redshift_series_plot
compare_simulations_plot
density_histogram_plot
density_profile_plot
metallicity_profile_plot
mass_profile_plot
sfr_txt_plot
temperature_histogram_plot
rho_temp_plot
kennicutt_schmidt_plot
```

### Testing

```@docs
comparison                    
deep_comparison
```

## Internal API

### Data Acquisition

```@docs
GADGETPlotting.get_snapshot_path
GADGETPlotting.get_time_evolution
GADGETPlotting.get_position
GADGETPlotting.get_density
GADGETPlotting.get_hsml
GADGETPlotting.get_mass
GADGETPlotting.get_metallicity
GADGETPlotting.get_temperature
GADGETPlotting.get_age
GADGETPlotting.get_birth_place
GADGETPlotting.get_sfr_txt
```

### Auxiliary

```@docs
GADGETPlotting.relative
GADGETPlotting.make_video
GADGETPlotting.smooth_window
GADGETPlotting.density_profile
GADGETPlotting.metallicity_profile
GADGETPlotting.mass_profile
GADGETPlotting.compute_cmdf
GADGETPlotting.kennicutt_schmidt_law
GADGETPlotting.format_error
GADGETPlotting.pass_all
GADGETPlotting.energy_integrand
GADGETPlotting.num_integrate
GADGETPlotting.center_of_mass
```

## Index

```@index
```