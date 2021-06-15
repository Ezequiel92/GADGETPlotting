using Base: func_for_method_checked, Float64
############################################################################################
# Auxiliary functions
############################################################################################

"""
    relative(
        p::Plots.Plot,
        rx::Float64,
        ry::Float64,
        rz::Union{Float64, Nothing} = nothing; 
        <keyword arguments>
    )::Union{NTuple{2, Float64}, NTuple{3, Float64}}
    
Give the absolute coordinates within a plot, from the relative ones.

If any of the axes are in a logarithmic scale, you should set `log` appropriately to get 
correct results.

# Arguments 
- `p::Plots.Plot`: Plot for which the absolute coordinates will be calculated.
- `rx::Float64`: relative x coordinate, `rx` ∈ [0, 1].
- `ry::Float64`: relative y coordinate, `ry` ∈ [0, 1].
- `rz::Union{Float64,Nothing} = nothing`: relative z coordinate, `rz` ∈ [0, 1].
- `log::Union{NTuple{2, Bool}, NTuple{3, Bool}} = (false, false, false)`: If the x, y or z 
  axes are in a logarithmic scale.

# Returns
- A Tuple with the absolute coordinates: (x, y) or (x, y, z).

# Examples
```julia-repl
julia> GADGETPlotting.relative(plot(rand(100)), 0.5, 0.5)
(50.5, 0.5047114800322484)

julia> GADGETPlotting.relative(surface(rand(100, 100)), 0.5, 0.5, 0.5)
(50.5, 50.5, 0.5000284432744244)
```
"""
function relative(
    p::Plots.Plot,
    rx::Float64,
    ry::Float64,
    rz::Union{Float64, Nothing} = nothing;
    log::Union{NTuple{2, Bool}, NTuple{3, Bool}} = (false, false, false),
)::Union{NTuple{2, Float64}, NTuple{3, Float64}}

    if log[1]
        xlims = log10.(Plots.xlims(p))
        ax = 10.0^(xlims[1] + rx * (xlims[2] - xlims[1]))
    else
        xlims = Plots.xlims(p)
        ax = xlims[1] + rx * (xlims[2] - xlims[1])
    end

    if log[2]
        ylims = log10.(Plots.ylims(p))
        ay = 10.0^(ylims[1] + ry * (ylims[2] - ylims[1]))
    else
        ylims = Plots.ylims(p)
        ay = ylims[1] + ry * (ylims[2] - ylims[1])
    end

    if rz === nothing

        return ax, ay

    else

        (
            length(log) == 3 || 
            error("If you have 3D coordinates, log has to have three values.")
        )

        if log[3]
            zlims = log10.(Plots.zlims(p))
            az = 10.0^(zlims[1] + rz * (zlims[2] - zlims[1]))
        else
            zlims = Plots.zlims(p)
            az = zlims[1] + rz * (zlims[2] - zlims[1])
        end

        return ax, ay, az

    end
end

"""
    make_video(
        source_path::String,
        source_format::String,
        output_path::String,
        output_filename::String,
        frame_rate::Int64,
    )::Nothing
	
Make an MP4 video from a series of images. 

The H.264 codec with no compression is used and the source images can be in
any format available in [ImageIO.jl](https://github.com/JuliaIO/ImageIO.jl).

# Arguments
- `source_path::String`: Path to the directory containing the images.	
- `source_format::String`: File format of the source images, e.g. ".png", ".svg", etc.
- `output_path::String`: Path to the directory where the resulting video will be saved.
- `output_filename::String`: Name of the video to be generated without extension.	
- `frame_rate::Int64`: Frame rate of the video to be generated.
"""
function make_video(
    source_path::String,
    source_format::String,
    output_path::String,
    output_filename::String,
    frame_rate::Int64,
)::Nothing

    # Loads the target images
    imagestack = [load(image) for image in glob("*" * source_format, source_path)]

    (
        !isempty(imagestack) ||
        error("I couldn't find any '$source_format' images in '$source_path'.")
    )

    # Creates the video with the specified frame rate and filename
    VideoIO.save(
        joinpath(output_path, output_filename * ".mp4"),
        imagestack,
        framerate = frame_rate;
        encoder_options = (crf = 0, preset = "ultrafast"),
        codec_name = "libx264rgb",
    )

    return nothing
end

"""
    smooth_window(
        x_data::Vector{<:Real},
        y_data::Vector{<:Real},
        bins::Int64,
    )::NTuple{2, Vector{Float64}}

Separate the values in `x_data` in `bins` contiguous windows, and replaces every x and y 
value within the window with the corresponding mean to smooth out the data. 

# Arguments
- `x_data::Vector{<:Real}`: x-axis data.
- `y_data::Vector{<:Real}`: y-axis data.
- `bins::Int64`: Number of windows to be used in the smoothing.
- `log::Bool = false`: If the x-axis data will be separated using logarithmic bins.

# Returns
- A Tuple with two arrays containing the smoothed out x and y data.
"""
function smooth_window(
    x_data::Vector{<:Real},
    y_data::Vector{<:Real},
    bins::Int64;
    log::Bool = false,
)::NTuple{2, Vector{Float64}}

    # Dimension consistency check
    (
        length(x_data) == length(y_data) ||
        throw(DimensionMismatch("The input vectors should have the same length."))
    )

    if log 
        # First positive value of the x axis
        start = log10(minimum(x -> x <= 0.0 ? Inf : x, x_data))
        # Logarithmic widths of the smoothing windows
        width = (log10(maximum(x_data)) - start) / bins
    else
        # First value of the x axis
        start = minimum(x_data)
        # Linear widths of the smoothing windows
        width = (maximum(x_data) - start) / bins
    end

    # Initialize output arrays.
    smooth_x_data = Vector{Float64}(undef, bins)
    smooth_y_data = Vector{Float64}(undef, bins)

    @inbounds for i in eachindex(smooth_x_data, smooth_y_data)

        # Find the indices of `x_data` which fall within window `i`
        if log
            idx = findall(
                x -> 10.0^(start + width * (i - 1)) <= x < 10.0^(start + width * i), 
                x_data,
            )
        else 
            idx = findall(x -> start + width * (i - 1) <= x < start + width * i, x_data)
        end
		
		if isempty(idx)
			error("Using $bins bins is too high for the data, lower it.")
		else
			# Store mean values in output arrays
			smooth_x_data[i] = sum(x_data[idx]) / length(idx)
			smooth_y_data[i] = sum(y_data[idx]) / length(idx)
		end
    end

    return smooth_x_data, smooth_y_data
end

@doc raw"""
    density_profile(
        mass_data::Vector{Float64},
        distance_data::Vector{Float64},
        max_radius::Float64,
        bins::Int64,
    )::NTuple{2, Vector{Float64}}
	
Compute a density profile up to a radius `max_radius`. 

It divides a sphere of radius `max_radius`, centered at (0, 0, 0), in `bins` spherical 
shells of equal width `max_radius` / `bins`. This results in a volume for the ``n``-th shell of

```math
V_n = \frac{4}{3}\,\pi\,\mathrm{width}^3\,(3\,n^2 - 3\,n + 1) \, .
```

So, the assigned density for that shell is ``\rho_n = M_n / V_n`` where ``M_n`` is the total 
mass within the shell.

`max_radius` and `distance_data` must be in the same length units.

# Arguments
- `mass_data::Vector{Float64}`: Masses of the particles.
- `distance_data::Vector{Float64}`: Radial distances of the particles. 
- `max_radius::Float64`: Maximum distance up to which the profile will be calculated.
- `bins::Int64`: Number of subdivisions of [0, `max_radius`] to be used for the profile.

# Returns
- A Tuple with two arrays. 
  The first with the radial distances and the second with the densities, of each shell.
"""
function density_profile(
    mass_data::Vector{Float64},
    distance_data::Vector{Float64},
    max_radius::Float64,
    bins::Int64,
)::NTuple{2, Vector{Float64}}

    # Dimension consistency check
    (
        length(mass_data) == length(distance_data) ||
        throw(DimensionMismatch("The input vectors should have the same length."))
    )

    # Width of each spherical shell used to calculate the density
    width = max_radius / bins

    # Initialize output arrays
    x_data = Vector{Float64}(undef, bins)
    y_data = Vector{Float64}(undef, bins)

    @inbounds for i in eachindex(x_data, y_data)
        # Find the indices of `distance_data` which fall within window `i`
        idx = findall(x -> width * (i - 1) <= x < width * i, distance_data)

        if isempty(idx)
            x_data[i] = width * (i - 0.5)
            y_data[i] = 0.0
        else
            total_mass = sum(mass_data[idx])
            volume = 4.0 / 3.0 * π * width^3.0 * (3.0 * i * i - 3.0 * i + 1.0)

            # Mean distance for window i
            x_data[i] = sum(distance_data[idx]) / length(idx)
            # Density for window i
            y_data[i] = total_mass / volume
        end
    end

    return x_data, y_data
end

@doc raw"""
    metallicity_profile(
        mass_data::Vector{Float64},
        distance_data::Vector{Float64},
        z_data::Vector{Float64},
        max_radius::Float64,
        bins::Int64,
    )::NTuple{2, Vector{Float64}}
	
Compute a metallicity profile up to a radius `max_radius` (normalized to solar metallicity).

It divides a sphere of radius `max_radius`, centered at (0, 0, 0), in `bins` spherical 
shells of equal width `max_radius / bins`. This results in a relative metallicity for 
the ``n``-th shell of

```math
\rho_n = \dfrac{z_n}{M_n\,Z_\odot} \, ,
```

where ``M_n`` is the total mass and ``z_n`` the total content of metals, within the shell.

`max_radius` and `distance_data` must be in the same length units, 
and `z_data` and `mass_data` must be in the same mass units.

# Arguments
- `mass_data::Vector{Float64}`: Masses of the particles.
- `distance_data::Vector{Float64}`: Radial distances of the particles. 
- `z_data::Vector{Float64}`: Metal content of the particles in mass units.
- `max_radius::Float64`: Maximum distance up to which the profile will be calculated.
- `bins::Int64`: Number of subdivisions of [0, `max_radius`] to be used for the profile.

# Returns
- A Tuple of two arrays.
  The first with the radial distances and the second with the metallicities, of each shell.
"""
function metallicity_profile(
    mass_data::Vector{Float64},
    distance_data::Vector{Float64},
    z_data::Vector{Float64},
    max_radius::Float64,
    bins::Int64,
)::NTuple{2, Vector{Float64}}

    # Dimension consistency check
    (
        length(mass_data) == length(distance_data) == length(z_data) ||
        throw(DimensionMismatch("The input vectors should have the same length."))
    )

    # Width of each spherical shell used to calculate the metallicity
    width = max_radius / bins

    # Initialize output arrays
    x_data = Vector{Float64}(undef, bins)
    y_data = Vector{Float64}(undef, bins)

    @inbounds for i in eachindex(x_data, y_data)
        # Find the indices of `distance_data` which fall within window `i`
        idx = findall(x -> width * (i - 1) <= x < width * i, distance_data)

        if isempty(idx)
            x_data[i] = width * (i - 0.5)
            y_data[i] = 0.0
        else
            total_mass = sum(mass_data[idx])
            total_z = sum(z_data[idx])

            # Mean distance for window i
            x_data[i] = sum(distance_data[idx]) / length(idx)
            # Metallicity for window i
            y_data[i] = (total_z / total_mass) / SOLAR_METALLICITY
        end
    end

    return x_data, y_data
end

@doc raw"""
    mass_profile(
        mass_data::Vector{Float64},
        distance_data::Vector{Float64},
        max_radius::Float64,
        bins::Int64,
    )::NTuple{2, Vector{Float64}}
	
Compute an cumulative mass profile up to a radius `max_radius`. 

It divides a sphere of radius `max_radius`, centered at (0, 0, 0), in `bins` spherical 
shells of equal width `max_radius / bins`. This results in a cumulative mass for 
the ``n``-th shell of

```math
M_n = \sum_{i = 1}^n m_i \, ,
```

where ``m_i`` is the total mass within the ``i``-th shell.

`max_radius` and `distance_data` must be in the same length units.

# Arguments
- `mass_data::Vector{Float64}`: Masses of the particles.
- `distance_data::Vector{Float64}`: Radial distances of the particles. 
- `max_radius::Float64`: Maximum distance up to which the profile will be calculated.
- `bins::Int64`: Number of subdivisions of [0, `max_radius`] to be used for the profile.

# Returns
- A Tuple of two arrays.
  The first with the radial distances and the second with the accumulated masses, 
  of each shell.
"""
function mass_profile(
    mass_data::Vector{Float64},
    distance_data::Vector{Float64},
    max_radius::Float64,
    bins::Int64,
)::NTuple{2, Vector{Float64}}

    # Dimension consistency check
    (
        length(mass_data) == length(distance_data) ||
        throw(DimensionMismatch("The input vectors should have the same length."))
    )

    # Width of each spherical shell used to calculate the mass
    width = max_radius / bins

    # Initialize output arrays
    x_data = Vector{Float64}(undef, bins)
    y_data = Vector{Float64}(undef, bins)

    @inbounds for i in eachindex(x_data, y_data)
        # Indices of `distance_data` within window i
        idx = findall(x -> width * (i - 1) <= x < width * i, distance_data)
		
		if isempty(idx)
            x_data[i] = width * (i - 0.5)
        else
			# Mean distance for window i
			x_data[i] = sum(distance_data[idx]) / length(idx)
		end
		
		# Mass for window i
		y_data[i] = sum(mass_data[idx])
    end

    return x_data, cumsum(y_data)
end

@doc raw"""
    compute_cmdf(
        mass_data::Vector{Float64},
        metallicity_data::Vector{Float64},
        max_Z::Float64,
        bins::Int64; 
        <keyword arguments>
    )::NTuple{2, Vector{Float64}}
	
Compute the cumulative metallicity distribution function up to a metallicity of `max_Z`. 

The CMDF is calculated separating in `bins` windows the stellar metallicity, within the 
range [0, `max_Z`]. For the ``n``-th window the CMDF is

```math
\sum_{i = 1}^n \dfrac{m_i}{M_T} \quad \mathrm{vs.} \quad \bar{Z}_n \, ,
```

or, for `x_norm = true`, 

```math
\sum_{i = 1}^n \dfrac{m_i}{M_T} \quad \mathrm{vs.} \quad \dfrac{\bar{Z}_n}{\mathrm{max}(\bar{Z}_n)} \, ,
```

where ``M_T`` is the total stellar mass, ``m_i`` the stellar mass of the ``i``-th window and
``\bar{Z}_n`` the mean stellar metallicity of the ``n``-th window.

`mass_data` and `metallicity_data` must be in the same mass units.

# Arguments
- `mass_data::Vector{Float64}`: Masses of the particles.
- `metallicity_data::Vector{Float64}`: Metallicities of the particles. 
- `max_Z::Float64`: Maximum metallicity up to which the CMDF will be calculated.
- `bins::Int64`: Number of subdivisions of [0, `max_Z`] to be used for the CMDF.
- `x_norm::Bool = false`: If the x axis will be normalized to its maximum value. 

# Returns
- A Tuple of two Arrays.
  The first with the metallicities and the second with the accumulated masses, 
  of each window.
"""
function compute_cmdf(
    mass_data::Vector{Float64},
    metallicity_data::Vector{Float64},
    max_Z::Float64,
    bins::Int64;
    x_norm::Bool = false,
)::NTuple{2, Vector{Float64}}

    # Dimension consistency check
    (
        length(mass_data) == length(metallicity_data) ||
        throw(DimensionMismatch("The input vectors should have the same length."))
    )
	
	# Dimensionless metallicity
	Z = metallicity_data ./ mass_data

    # If required, normalize the x axis
    if x_norm
        Z = Z ./ max_Z
        # Width of the metallicity bins
        width = 1.0 / bins
    else
        # Width of the metallicity bins
	    width = max_Z / bins  
    end
	
	# Total star mass
	total_m = sum(mass_data)
	
	# Initialize output arrays
    x_data = Vector{Float64}(undef, bins)
    y_data = Vector{Float64}(undef, bins)

    @inbounds for i in eachindex(x_data, y_data)
		idx = findall(x -> width * (i - 1) <= x < width * i, Z)
		
		if isempty(idx)
            x_data[i] = width * (i - 0.5)
        else
			# Mean metallicity for window i
			x_data[i] = sum(Z[idx]) / length(idx)
		end
		
		# Mass fraction for window i
		y_data[i] = sum(mass_data[idx]) / total_m			
		
	end
	
	return x_data, cumsum(y_data)
end

@doc raw"""
    kennicutt_schmidt_law(
        gas_mass_data::Vector{Float64},
        gas_distance_data::Vector{Float64},
        temperature_data::Vector{Float64},
        star_mass_data::Vector{Float64},
        star_distance_data::Vector{Float64},
        age_data::Vector{Float64},
        temp_filter::Float64,
        age_filter::Float64,
        max_r::Float64; 
        <keyword arguments>
    )::Union{Nothing, Dict{String, Any}}
	
Compute the mass area density and the SFR area density for the [Kennicutt-Schmidt law](https://en.wikipedia.org/wiki/Kennicutt%E2%80%93Schmidt_law). 

The area densities are calculated by projecting the positions of the stars and of the 
particles of gas to the x-y plane. Then the space is subdivided in `bins` concentric 
rings from 0 to `max_r`, each of equal width `max_r` / `bins`. This results in an area 
for the ``n``-th ring of

```math
A_n = π \, \mathrm{width}^2 \, (2 \, n - 1) \, .
```

So, the assigned SFR for that ring is 

```math
\Sigma_\mathrm{SFR}^n = \frac{M_*^n}{\mathrm{age\_filter}\,A_n} \, ,
```

where ``M_*^n`` is the total mass of stars younger than `age_filter` within the ring.

Equivalently, the mass area density of the gas is given by

```math
\Sigma_\rho^n = \frac{M_\rho^n}{A_n} \, ,
```

where ``M_\rho^n`` is the total mass of gas colder than `temp_filter` within the ring.

`temp_filter` and `temperature_data` must be in the same temperature units, 
and `age_filter` and `age_data` must be in the same time units.

# Arguments
- `gas_mass_data::Vector{Float64}`: Masses of the gas particles.
- `gas_distance_data::Vector{Float64}`: 2D radial distances of the gas particles. 
- `temperature_data::Vector{Float64}`: Temperatures of the gas particles.
- `star_mass_data::Vector{Float64}`: Masses of the stars.
- `star_distance_data::Vector{Float64}`: 2D radial distances of the stars.
- `age_data::Vector{Float64}`: Ages of the stars.
- `temp_filter::Float64`: Maximum temperature allowed for the gas particles.
- `age_filter::Unitful.Quantity`: Maximum age allowed for the stars.
- `max_r::Float64`: Maximum distance up to which the parameters will be calculated.
- `bins::Int64 = 50`: Number of subdivisions of [0, `max_r`] to be used. It has to be at 
  least 5.

# Returns
- A dictionary with three entries.
  - Key "RHO" => Logarithm of the area mass densities.
  - Key "SFR" => Logarithm of the SFR area densities.
  - Key "LM" => Linear model given by [GLM.jl](https://github.com/JuliaStats/GLM.jl).
- Or `nothing` if there are less than 5 data points in the end result.
"""
function kennicutt_schmidt_law(
    gas_mass_data::Vector{Float64},
    gas_distance_data::Vector{Float64},
    temperature_data::Vector{Float64},
    star_mass_data::Vector{Float64},
    star_distance_data::Vector{Float64},
    age_data::Vector{Float64},
    temp_filter::Float64,
    age_filter::Float64,
    max_r::Float64;
    bins::Int64 = 50,
)::Union{Nothing, Dict{String, Any}}

    # Bin size check
    if bins < 5
        error("You have to use at least 5 bins.")
    end

    # Filter out hot gas particles
    cold_gas_mass = deleteat!(copy(gas_mass_data), temperature_data .> temp_filter)
    cold_gas_distance = deleteat!(copy(gas_distance_data), temperature_data .> temp_filter)

    # Filter out old stars
    young_star_mass = deleteat!(copy(star_mass_data), age_data .> age_filter)
    young_star_distance = deleteat!(copy(star_distance_data), age_data .> age_filter)

    r_width = max_r / bins

    # Initialize output arrays
    x_data = Vector{Float64}(undef, bins)
    y_data = Vector{Float64}(undef, bins)

    @inbounds for i in eachindex(x_data, y_data)

        # Gas.
		idx_gas = findall(x ->  r_width * (i - 1) <= x < r_width * i, cold_gas_distance)
        gas_mass = sum(cold_gas_mass[idx_gas])
		# Gas area density for window i
		x_data[i] = gas_mass / (π * r_width * r_width * (2.0 * i - 1.0))

        # Stars.
        idx_star = findall(x ->  r_width * (i - 1) <= x < r_width * i, young_star_distance)
        sfr = sum(young_star_mass[idx_star]) / age_filter 
        # SFR area density for window i
        y_data[i] = sfr / (π * r_width * r_width * (2.0 * i - 1.0))		
		
	end

    # Filter out zeros
    deleteat!(y_data, x_data .<= 0.0)
    filter!(x -> x > 0.0, x_data)
    deleteat!(x_data, y_data .<= 0.0)
    filter!(y -> y > 0.0, y_data)

    # Set logarithmic scaling
    x_data = log10.(x_data)
    y_data = log10.(y_data)

    # If there are less than 5 data points return nothing
    if length(x_data) < 5
        return nothing
    end

    # Compute linear fit
    X = [ones(length(x_data)) x_data]
    linear_model = lm(X, y_data)

    return Dict("RHO" => x_data, "SFR" => y_data, "LM" => linear_model)
end

@doc raw"""
    quantities_2D(
        gas_mass_data::Vector{Float64},
        gas_distance_data::Vector{Float64},
        temperature_data::Vector{Float64},
        star_mass_data::Vector{Float64},
        star_distance_data::Vector{Float64},
        age_data::Vector{Float64},
        metal_mass_data::Matrix{Float64},
        fmol_data::Vector{Float64},
        temp_filter::Float64,
        age_filter::Float64,	
        max_r::Float64;
        bins::Int64 = 50,
    )::Union{Nothing, Dict{String, Any}}
	
Compute the surface density of several quantities. 

The area densities are calculated by projecting the positions of the stars and the 
particles of gas to the x-y plane. Then the space is subdivided in `bins` concentric 
rings from 0 to `max_r`, each of equal width `max_r` / `bins`. This results in an area 
for the ``n``-th ring of

```math
A_n = π \, \mathrm{width}^2 \, (2 \, n - 1) \, .
```

So, the assigned SFR for that ring is 

```math
\Sigma_\mathrm{SFR}^n = \frac{M_*^n}{\mathrm{age\_filter}\,A_n} \, ,
```

where ``M_*^n`` is the total mass of stars younger than `age_filter` within the ring.

Equivalently, the mass area density of the gas and stars is given by

```math
\Sigma_\rho^n = \frac{M_\rho^n}{A_n} \, ,
```
```math
\Sigma_*^n = \frac{M_*^n}{A_n} \, ,
```

where ``M_\rho^n`` is the total mass of gas colder than `temp_filter` within the ring.

`temp_filter` and `temperature_data` must be in the same temperature units, 
and `age_filter` and `age_data` must be in the same time units.

The rest of the parameters are define as follows

```math
\mathrm{SSFR}^n = \frac{\Sigma_\mathrm{SFR}^n}{\Sigma_{*\mathrm{(total)}}^n} \, ,
```
```math
\mathrm{SFE}^n = \frac{\Sigma_\mathrm{SFR}^n}{\Sigma_{\rho\mathrm{(total)}}^n} \, ,
```
```math
\mathrm{P}^n = \Sigma_{\rho\mathrm{(total)}}^n\left(\Sigma_{\rho\mathrm{(total)}}^n + \Sigma_{*\mathrm{(total)}}^n\right) \, ,
```
```math
\mathrm{Ψ/H_2}^n = \frac{\Sigma_\mathrm{SFR}^n * A_n}{M_{H_2}^n} \, ,
```

where ``\Sigma_{\rho\mathrm{(total)}}^n`` is the mass density of all the gas, 
``\Sigma_{*\mathrm{(total)}}^n`` is the stellar density of all the stars and
``M_{H_2}^n`` is the total mass of molecular Hydrogen.

# Arguments
- `gas_mass_data::Vector{Float64}`: Masses of the gas particles.
- `gas_distance_data::Vector{Float64}`: 2D radial distances of the gas particles. 
- `temperature_data::Vector{Float64}`: Temperatures of the gas particles.
- `star_mass_data::Vector{Float64}`: Masses of the stars.
- `star_distance_data::Vector{Float64}`: 2D radial distances of the stars.
- `age_data::Vector{Float64}`: Ages of the stars.
- `metal_mass_data::Matrix{Float64}`: Masses of the individual elements (H, O, He, C, etc.)
  within the gas particles.
- `fmol_data::Vector{Float64}`: Fraction of molecular Hydrogen of the gas particles.
- `temp_filter::Float64`: Maximum temperature allowed for the gas particles.
- `age_filter::Unitful.Quantity`: Maximum age allowed for the stars.
- `max_r::Float64`: Maximum distance up to which the parameters will be calculated.
- `bins::Int64 = 50`: Number of subdivisions of [0, `max_r`] to be used. It has to be at 
  least 5.

# Returns
- A dictionary with nine entries.
  - Key "GAS" => Surface mass density of gas
  - Key "COLD_GAS" => Surface mass density of cold gas
  - Key "STARS" => Surface mass density of stars
  - Key "OH" => 12 + log10(oxygen_mass / hydrogen_mass)
  - Key "SFR" => Star formation rate surface density
  - Key "SSFR" => Specific star formation rate surface density
  - Key "SFE" => Star formation efficiency surface density
  - Key "P" => Proportional to the pressure
  - Key "Ψ_FMOL" => Star formation rate per unit of molecular gas
"""
function quantities_2D(
    gas_mass_data::Vector{Float64},
    gas_distance_data::Vector{Float64},
    temperature_data::Vector{Float64},
    star_mass_data::Vector{Float64},
    star_distance_data::Vector{Float64},
    age_data::Vector{Float64},
	metal_mass_data::Matrix{Float64},
    fmol_data::Vector{Float64},
    temp_filter::Float64,
    age_filter::Float64,	
    max_r::Float64;
    bins::Int64 = 50,
)::Union{Nothing, Dict{String, Any}}

    # Bin size check
    if bins < 5
        error("You have to use at least 5 bins.")
    end
	
	# Oxygen
	O = metal_mass_data[4, :]
	
	# Hydrogen
	H = metal_mass_data[7, :]

    r_width = max_r / bins

    # Initialize output arrays
	GAS = Vector{Float64}(undef, bins)
    COLD_GAS = Vector{Float64}(undef, bins)
	STARS = Vector{Float64}(undef, bins)
	OH = Vector{Union{Float64, Missing}}(undef, bins)
    SFR = Vector{Float64}(undef, bins)
    SSFR = Vector{Union{Float64, Missing}}(undef, bins)
	SFE = Vector{Union{Float64, Missing}}(undef, bins)
	P = Vector{Float64}(undef, bins)
    Ψ_FMOL = Vector{Union{Float64, Missing}}(undef, bins)

    @inbounds for i in 1:bins

        # Gas
		
		# Indices of gas particles within the window i
		idx_gas = findall(x ->  r_width * (i - 1) <= x < r_width * i, gas_distance_data)
		# Masses of gas particles within the window i
		gas_masses = gas_mass_data[idx_gas]
		# Total gas mass within the window i
		gas_mass = sum(gas_masses)	
		# Total cold gas mass within the window i
		cold_gas_mass = sum(deleteat!(gas_masses, temperature_data[idx_gas] .> temp_filter))
		
		# Metals 
		
		# Total oxygen mass
		m_O = sum(O[idx_gas])
		# Total Hydrogen mass
		m_H = sum(H[idx_gas])

        # Molecular Hydrogen

        m_H_mol = sum(H[idx_gas] .* fmol_data[idx_gas])
		
		# Stars
		
		# Indices of stellar particles within the window i
		idx_star = findall(x ->  r_width * (i - 1) <= x < r_width * i, star_distance_data)
		# Masses of stellar particles within the window i
		star_masses = star_mass_data[idx_star]
		# Total stellar mass within the window i
		star_mass = sum(star_masses)	
		# Total young stellar mass within the window i
		young_star_mass = sum(deleteat!(star_masses, age_data[idx_star] .> age_filter))
		# Star formation rate for the window i
		sfr = young_star_mass / age_filter 
		
	    # Ring area
	    a = π * r_width * r_width * (2.0 * i - 1.0)
		
		# Quantities
			
		GAS[i] = gas_mass / a                  # Surface mass density of gas
        COLD_GAS[i] = cold_gas_mass / a        # Surface mass density of cold gas
		STARS[i] = star_mass / a               # Surface mass density of stars
		SFR[i] = sfr / a                       # Star formation rate surface density
        P[i] = GAS[i] * (GAS[i] + STARS[i])    # Proportional to the pressure

        # Specific star formation rate surface density
        if STARS[i] == 0.0                     
            SSFR[i] = missing
        else
		    SSFR[i] = SFR[i] / STARS[i]           
        end

        # Star formation efficiency surface density
        if GAS[i] == 0.0                       
            SFE[i] = missing
        else
		    SFE[i] = SFR[i] / GAS[i]           
        end

        # Star formation rate per unit of molecular gas
        if m_H_mol == 0.0                      
            Ψ_FMOL[i] = missing
        else
		    Ψ_FMOL[i] = sfr / m_H_mol           
        end 
         
        # Oxygen-Hydrogen relation = 12 + log10(oxygen_mass / hydrogen_mass)             
		if m_O == 0.0 || m_H == 0.0            
			OH[i] = missing
		else
			OH[i] = 12 + log10(m_O / m_H)
		end
		
	end

    return Dict(
		"GAS"      => GAS,	
        "COLD_GAS" => COLD_GAS,
		"STARS"    => STARS,
		"OH"       => OH,
		"SFR"      => SFR,
		"SSFR"     => SSFR,
		"SFE"      => SFE,
		"P"        => P,
        "Ψ_FMOL"   => Ψ_FMOL,
	)
end

"""
    format_error(mean::Float64, error::Float64)::String

Format the mean and error values.

It follows the traditional rules for error presentation. The error has only one significant  
digit, unless such digit is a one, in which case, two significant digits are used.  
The mean will have a number of digits such as to match the last significant position 
of the error. 

# Arguments 
- `mean::Float64`: Mean value.
- `error::Float64`: Error value. It must be positive.

# Returns
- A Tuple with the formatted mean and error values.

# Examples
```julia-repl
julia> format_error(69.42069, 0.038796)
(69.42, 0.04)

julia> format_error(69.42069, 0.018796)
(69.421, 0.019)

julia> format_error(69.42069, 0.0)
(69.42069, 0.0)

julia> format_error(69.42069, 73.4)
(0.0, 70.0)
```
"""
function format_error(mean::Float64, error::Float64)::NTuple{2, Float64}

    # Positive error check
    error >= 0.0 || error("The error must be a positive number.")

    if error == 0.0
        round_mean = mean
        round_error = error
    else
        sigdigit_pos = abs(log10(abs(error)))

        if error < 1.0
            if abs(mean) < error
                extra = 0
                round_mean = 0.0
            else
                first_digit = trunc(error * 10.0^(floor(sigdigit_pos) + 1.0))
                first_digit == 1.0 ? extra = 1 : extra = 0

                digits = ceil(Int64, sigdigit_pos) + extra
                round_mean = round(mean; digits)
            end
        else
            if abs(mean) < error
                extra = 0
                round_mean = 0.0
            else
                first_digit = trunc(error / 10.0^(floor(sigdigit_pos)))
                first_digit == 1.0 ? extra = 2 : extra = 1

                sigdigits = ceil(Int64, log10(abs(mean))) - ceil(Int64, sigdigit_pos) + extra
                round_mean = round(mean; sigdigits)
            end
        end

        round_error = round(error, sigdigits = 1 + extra)
    end

    return round_mean, round_error
end

"""
    pass_all(snap_file::String, type::String)::Vector{Int64}

Default filter function for the [read\\_blocks\\_over\\_all\\_files](https://ludwigboess.github.io/GadgetIO.jl/stable/api/#GadgetIO.read_blocks_over_all_files-Tuple{String,%20Array{String,%20N}%20where%20N}) function.

It does not filter out any particles, allowing the data acquisition functions to gather 
all the data. 

# Arguments 
- `snap_file::String`: Snapshot file path.
- `type::String`: Particle type.
    * "gas" -> Gas particle. 
    * "dark_matter" -> Dark matter particle.
    * "stars" -> Star particle.

# Returns
- A Vector with the indices of the allowed particles.
"""
function pass_all(snap_file::String, type::String)::Vector{Int64}

    # Select type of particle
    if type == "gas"
        type_num = 1
    elseif type == "dark_matter"
        type_num = 2
    elseif type == "stars"
        type_num = 5
    else
        error("Particle type '$type' not supported. 
        The supported types are 'gas', 'dark_matter' and 'stars'")
    end

    header = read_header(snap_file)

    return collect(1:header.npart[type_num])
end

@doc raw"""
    energy_integrand(header::GadgetIO.SnapshotHeader, a::Float64)::Float64

Give the integrand of the scale factor to physical time function

```math
\frac{1}{H\,\sqrt{\epsilon}} \, ,
``` 

where 

```math
\epsilon = \Omega_\lambda + \frac{1 - \Omega_\lambda - \Omega_0}{a^2} + \frac{\Omega_0}{a^3} \, , 
```
```math
H = H_0 \, a \, .
```

# Arguments 
- `header::GadgetIO.SnapshotHeader`: Header of the relevant snapshot file.
- `a::Float64`: Dimensionless scale factor.

# Returns
- The integrand in Gyr evaluated in `a` .
"""
function energy_integrand(header::GadgetIO.SnapshotHeader, a::Float64)::Float64

    # Ω_K (curvature)
    omega_K = 1.0 - header.omega_0 - header.omega_l
    # Energy function
    E = header.omega_0 / (a * a * a) + omega_K / (a * a) + header.omega_l
    # Hubble constant in 1 / Gyr
    H = header.h0 * HUBBLE_CONST * a

    # Integrand of the time integral in Gyr
    return 1.0 / (H * sqrt(E))
end

@doc raw"""
    num_integrate(
        func::Function, 
        inf_lim::Float64, 
        sup_lim::Float64, 
        steps::Int64 = 200,
    )::Float64

Give the numerical integral of `func` between `inf_val` and `sup_val`. 

The result is given by

```math
\int_\mathrm{inf\_lim}^\mathrm{sup\_lim} f(x) \mathrm{dx} \approx \sum_{i = 1}^\mathrm{steps} f(\mathrm{inf\_lim} + \mathrm{width}\,i ) \, ,
```

where ``\mathrm{width} = (\mathrm{sup\_lim} - \mathrm{inf\_lim}) / \mathrm{steps}``.

# Arguments 
- `func::Function`: 1D function to be integrated.
- `inf_lim::Float64`: Lower limit of the integral.
- `sup_lim::Float64`: Upper limit if the integral.
- `steps::Int64`: Number of subdivisions to be used for the discretization of 
  the [`inf_lim`, `sup_lim`] region.

# Returns
- The value of the integral.

# Examples
```julia-repl
julia> num_integrate(sin, 0, 3π)
1.9996298761360816

julia> num_integrate(x -> x^3 + 6 * x^2 + 9 * x + 2, 0, 4.69)
438.9004836958452

julia> num_integrate(x -> exp(x^x), 0, 1.0)
2.1975912134624904
```
"""
function num_integrate(
    func::Function, 
    inf_lim::Real, 
    sup_lim::Real, 
    steps::Int64 = 200,
)::Float64
    
    # Width of a single subinterval
    width = (sup_lim - inf_lim) / steps
    # Integrand evaluated at the rightmost value of each subinterval
    integrand = func.([inf_lim + width * i for i in 1:steps])

    # Final result of the numerical integration
    return sum(width .* integrand)
end

@doc raw"""
    center_of_mass(
        position_data::Matrix{<:Real},
        mass_data::Vector{<:Real},
    )::Union{NTuple{3, Float64}, Nothing}

Compute the center of mass as

```math
R_c = \frac{1}{M} \sum_n m_n \, r_n \, ,
```

where ``M = \sum_n m_n`` and ``m_n`` and ``r_n`` are the mass and distance 
from the origin of the ``n``-th particle.

If the length of ``R`` is less than ``10^-3`` the length of the larger position vector in
`position_data`, nothing is returned.

# Arguments
- `position_data::Matrix{<:Real}`: Positions of the particles.
- `mass_data::Vector{<:Real}`: Masses of the particles.

# Returns
- The center of mass in the unis of `position_data`.
"""
function center_of_mass(
    position_data::Matrix{<:Real},
    mass_data::Vector{<:Real},
)::Union{NTuple{3, Float64}, Nothing}

    # Total mass
    M = sum(mass_data)

    R = [0.0, 0.0, 0.0]
    for (col, mass) in zip(eachcol(position_data), mass_data)
        R .+= col .* mass
    end

    center_of_mass = (R[1] / M, R[2] / M, R[3] / M)

    if norm(center_of_mass) < max_length(position_data) * 10^3
        return nothing
    else
        return center_of_mass
    end

end

"""
    max_length(data::Matrix{<:Real})::Float64

Maximum norm of the positions in `data`.

`data` must be a matrix with three rows (x, y, and z coordinates respectively) and where 
each column is a position.

# Arguments
- `data::Matrix{<:Real}`: Positions of the particles.

# Returns
- The maximum norm of the position vectors in `data`.
"""
@inline function max_length(data::Matrix{<:Real})::Float64

    return maximum([norm(col) for col in eachcol(data)])

end

"""
    comparison(
        x::Union{Real, AbstractArray{<:Real}, Tuple{Vararg{Real}}}, 
        y::Union{Real, AbstractArray{<:Real}, Tuple{Vararg{Real}}}; 
        atol::Float64 = 1e-5, 
        rtol::Float64 = 1e-5,
    )::Bool

Determine if two numbers, numeric arrays, or numeric tuples are approximately equal.

# Arguments
- `x::Union{Real, AbstractArray{<:Real}, Tuple{Vararg{Real}}}`: First element to be compared.
- `y::Union{Real, AbstractArray{<:Real}, Tuple{Vararg{Real}}}`: Second element to be compared.
- `atol::Float64 = 1e-5`: Absolute tolerance.
- `rtol::Float64 = 1e-5`: Relative tolerance.

# Returns
- Return `true` if every pair of elements (X, Y) in (`x`, `y`) pass
  ```julia
  norm(X - Y) <= max(atol, rtol * max(norm(X), norm(Y)))
  ````
"""
function comparison(
    x::Union{Real, AbstractArray{<:Real}, Tuple{Vararg{Real}}}, 
    y::Union{Real, AbstractArray{<:Real}, Tuple{Vararg{Real}}}; 
    atol::Float64 = 1e-5, 
    rtol::Float64 = 1e-5,
)::Bool

    return all(isapprox.(x, y; atol, rtol))

end

"""
    comparison(x, y; atol::Float64 = 1e-5, rtol::Float64 = 1e-5)::Bool

Determine if two elements are equal, as per the [isequal](https://docs.julialang.org/en/v1/base/base/#Base.isequal) function.

# Arguments
- `x`: First element to be compared.
- `y`: Second element to be compared.
- `atol::Float64 = 1e-5`: Absolute tolerance (for compatibility).
- `rtol::Float64 = 1e-5`: Relative tolerance (for compatibility).

# Returns
- Returns `isequal(x, y)`.
"""
function comparison(x, y; atol::Float64 = 1e-5, rtol::Float64 = 1e-5)::Bool

    return isequal(x, y)

end

"""
    deep_comparison(
        x::Dict, 
        y::Dict; 
        atol::Float64 = 1e-5, 
        rtol::Float64 = 1e-5,
    )::Bool

Determine if two dictionaries are approximately equal.

Numeric elements are compared with the [`comparison`](@ref) function, everything else with 
the [isequal](https://docs.julialang.org/en/v1/base/base/#Base.isequal) function.

# Arguments
- `x::Dict`: First dictionary to be compared.
- `y::Dict`: Second dictionary to be compared.
- `atol::Float64 = 1e-5`: Absolute tolerance for numeric elements.
- `rtol::Float64 = 1e-5`: Relative tolerance for numeric elements.

# Returns
- Return `true` if every pair of elements within the dictionaries pass the equality tests.
"""
function deep_comparison(
    x::Dict, 
    y::Dict; 
    atol::Float64 = 1e-5, 
    rtol::Float64 = 1e-5,
)::Bool

    if keys(x) != keys(y)
        return false
    end

    return all([comparison(x[key], y[key]; atol, rtol) for key in keys(x)])
    
end

"""
    deep_comparison(
        x::Union{AbstractArray, Tuple}, 
        y::Union{AbstractArray, Tuple}; 
        atol::Float64 = 1e-5, 
        rtol::Float64 = 1e-5,
    )::Bool

Determines if two arrays or tuples are approximately equal.

Numeric elements are compared with the [`comparison`](@ref) function, everything else with 
the [isequal](https://docs.julialang.org/en/v1/base/base/#Base.isequal) function.

# Arguments
- `x::Union{AbstractArray, Tuple}`: First array to be compared.
- `y::Union{AbstractArray, Tuple}`: Second array to be compared.
- `atol::Float64 = 1e-5`: Absolute tolerance for numeric elements.
- `rtol::Float64 = 1e-5`: Relative tolerance for numeric elements.

# Returns
- Return `true` if every pair of elements pass the equality tests.
"""
function deep_comparison(
    x::Union{AbstractArray, Tuple}, 
    y::Union{AbstractArray, Tuple}; 
    atol::Float64 = 1e-5, 
    rtol::Float64 = 1e-5,
)::Bool

    if length(x) != length(y)
        return false
    end

    return all([comparison(X, Y; atol, rtol) for (X, Y) in zip(x, y)])
    
end

"""
    set_vertical_flags(
        flags::Union{Tuple{Vector{<:Real}, Vector{<:AbstractString}}, Nothing}, 
        plot::Plots.Plot; 
        <keyword arguments>
    )::Plots.Plot

Draw vertical lines at specified positions.

If you only want labels for the vertical lines, the original plot should have `label = ""`.

# Arguments
- `flags::Union{Tuple{Vector{<:Real}, Vector{<:AbstractString}}, Nothing} = nothing`: The first 
  vector in the Tuple has the positions of the vetical lines. The second has the 
  corresponding labels.
- `plot::Plots.Plot`: Plot to which the vertical lines will be added.

# Returns
- New plot with the vertical lines added.
"""
function set_vertical_flags(
    flags::Union{Tuple{Vector{<:Real}, Vector{<:AbstractString}}, Nothing}, 
    plot::Plots.Plot,
)::Plots.Plot

    if flags === nothing
        return plot
    end

    # Copy input data
    pl_out = deepcopy(plot)
    v_lines = copy(flags[1])
    line_labels = copy(flags[2])

    # Filter out values larger than the maximum of the original plot
	xlims = Plots.xlims(pl_out)
    deleteat!(line_labels, v_lines .> xlims[2])
    filter!(x -> x <= xlims[2], v_lines)

    if isempty(v_lines)
        return pl_out
    end
    
    # Draw vertical lines
    for (v_line, line_label) in zip(v_lines, line_labels)
		vline!(pl_out, [v_line], line = 2, label = line_label) 
	end

    plot!(
        pl_out,
        legend = :topright,
        legendfontsize = 22,
    )

    return pl_out
end