############################################################################################
# AUXILIARY FUNCTIONS.
############################################################################################

"""
    relative(
        p::Plots.Plot,
        rx::Float64,
        ry::Float64,
        rz::Union{Float64, Nothing} = nothing,
    )::Union{NTuple{2, Float64}, NTuple{3, Float64}}
    
Give the absolute coordinates for a Plot, given the relative ones.

# Arguments 
- `p::Plots.Plot`: Plot for which the absolute coordinates will be calculated.
- `rx::Float64`: relative x coordinate, rx ∈ [0, 1].
- `ry::Float64`: relative y coordinate, ry ∈ [0, 1].
- `rz::Union{Float64,Nothing} = nothing`: relative z coordinate, rz ∈ [0, 1].

# Returns
- A Tuple with the absolute coordinates: (x, y) or (x, y, z).
"""
function relative(
    p::Plots.Plot,
    rx::Float64,
    ry::Float64,
    rz::Union{Float64, Nothing} = nothing,
)::Union{NTuple{2, Float64}, NTuple{3, Float64}}

    # Plot axes limits.
    xlims = Plots.xlims(p)
    ylims = Plots.ylims(p)

    if rz === nothing
        return xlims[1] + rx * (xlims[2] - xlims[1]), ylims[1] + ry * (ylims[2] - ylims[1])
    else
        zlims = Plots.zlims(p)

        return xlims[1] + rx * (xlims[2] - xlims[1]),
        ylims[1] + ry * (ylims[2] - ylims[1]),
        zlims[1] + rz * (zlims[2] - zlims[1])
    end
end

"""
    makeVideo(
        source_path::String,
        source_format::String,
        output_path::String,
        output_filename::String,
        frame_rate::Int64,
    )::Nothing
	
Make a MP4 video from a series of images. 

The H.264 codec is used with no compression and the source images can be in
any format available in ImageIO.jl, e.g. ".png", ".svg", ".jpeg", etc.

# Arguments
- `source_path::String`: Path to the directory containing the images.	
- `source_format::String`: File format of the source images. 
- `output_path::String`: Path to the directory where the resulting video will be saved.
- `output_filename::String`: Name of the video to be generated without extension.	
- `frame_rate::Int64`: Frame rate of the video to be generated.
"""
function makeVideo(
    source_path::String,
    source_format::String,
    output_path::String,
    output_filename::String,
    frame_rate::Int64,
)::Nothing

    # Loads the target images.
    image_stack = [load(image) for image in glob("*" * source_format, source_path)]

    (
        !isempty(image_stack) ||
        error("I couldn't find any '$source_format' images in '$source_path'.")
    )

    # Creates the video with the specified frame rate and filename.
    properties = [:priv_data => ("crf" => "0", "preset" => "ultrafast")]
    encodevideo(
        output_path * output_filename * ".mp4",
        image_stack,
        framerate = frame_rate,
        AVCodecContextProperties = properties,
    )

    return nothing
end

"""
    smoothWindow(
        x_data::Vector{T} where {T <: Real},
        y_data::Vector{T} where {T <: Real},
        bins::Int64,
    )::NTuple{2, Vector{Float64}}

Separate the range of values of `x_data` in `bins` contiguous windows, and replaces 
every value within the window with the mean in order to smooth out the data. 

# Arguments
- `x_data::Vector{T} where {T <: Real}`: Data used to create the windows.
- `y_data::Vector{T} where {T <: Real}`: Data to be smoothed out.
- `bins::Int64`: Number of windows to be used in the smoothing.
- `log::Bool = false`: If the x axis will be divided using logarithmic bins.

# Returns
- The smooth data.
"""
function smoothWindow(
    x_data::Vector{T} where {T <: Real},
    y_data::Vector{T} where {T <: Real},
    bins::Int64;
    log::Bool = false,
)::NTuple{2, Vector{Float64}}

    # Dimension consistency check.
    (
        length(x_data) == length(y_data) ||
        throw(DimensionMismatch("The input vectors should have the same length."))
    )

    if log 
        # First positive value of the x axis.
        start = log10(minimum(x -> x <= 0 ? Inf : x, x_data))
        # Logarithmic widths of the smoothing windows.
        width = (log10(maximum(x_data)) - start) / bins
    else
        # First value of the x axis.
        start = minimum(x_data)
        # Linear widths of the smoothing windows.
        width = (maximum(x_data) - start) / bins
    end

    # Initialize output arrays.
    smooth_x_data = Vector{Float64}(undef, bins)
    smooth_y_data = Vector{Float64}(undef, bins)

    @inbounds for i in eachindex(smooth_x_data, smooth_y_data)

        # Find the indices of `x_data` which fall within window `i`.
        if log
            idx = findall(
                x -> 10^(start + width * (i - 1)) <= x < 10^(start + width * i), 
                x_data,
            )
        else 
            idx = findall(x -> start + width * (i - 1) <= x < start + width * i, x_data)
        end
		
		if isempty(idx)
			error("Using $bins bins is too high for the data, lower it.")
		else
			# Store mean values in output arrays.
			smooth_x_data[i] = sum(x_data[idx]) / length(idx)
			smooth_y_data[i] = sum(y_data[idx]) / length(idx)
		end
    end

    return smooth_x_data, smooth_y_data
end

"""
    densityProfile(
        mass_data::Vector{Float64},
        distance_data::Vector{Float64},
        max_radius::Float64,
        bins::Int64,
    )::NTuple{2, Vector{Float64}}
	
Compute a density profile up to a radius `max_radius`. 

`max_radius` and `distance_data` must be in the same units.

# Arguments
- `mass_data::Vector{Float64}`: Masses of the particles.
- `distance_data::Vector{Float64}`: Radial distances of the particles. 
- `max_radius::Float64`: Maximum distance up to which the profile will be calculated.
- `bins::Int64`: Number of subdivisions of [0, `max_radius`] to be used for the profile.

# Returns
- A Tuple of two Arrays. 
  The first with the radial distances and the second with the densities.
"""
function densityProfile(
    mass_data::Vector{Float64},
    distance_data::Vector{Float64},
    max_radius::Float64,
    bins::Int64,
)::NTuple{2, Vector{Float64}}

    # Dimension consistency check.
    (
        length(mass_data) == length(distance_data) ||
        throw(DimensionMismatch("The input vectors should have the same length."))
    )

    # Width of each spherical shell used to calculate the density.
    width = max_radius / bins

    # Initialize output arrays.
    x_data = Vector{Float64}(undef, bins)
    y_data = Vector{Float64}(undef, bins)

    @inbounds for i in eachindex(x_data, y_data)
        # Find the indices of `distance_data` which fall within window `i`.
        idx = findall(x -> width * (i - 1) <= x < width * i, distance_data)

        if isempty(idx)
            x_data[i] = width * (i - 0.5)
            y_data[i] = 0.0
        else
            total_mass = sum(mass_data[idx])
            volume = 4 / 3 * π * width^3 * (3 * i * i - 3 * i + 1)

            # Mean distance for window i.
            x_data[i] = sum(distance_data[idx]) / length(idx)
            # Density for window i.
            y_data[i] = total_mass / volume
        end
    end

    return x_data, y_data
end

"""
    metallicityProfile(
        mass_data::Vector{Float64},
        distance_data::Vector{Float64},
        z_data::Vector{Float64},
        max_radius::Float64,
        bins::Int64,
    )::NTuple{2, Vector{Float64}}
	
Compute a metallicity profile up to a radius `max_radius`, 
and normalize it to the solar metallicity.

`max_radius` and `distance_data` must be in the same units.
`z_data` and `mass_data` must be in the same units.

# Arguments
- `mass_data::Vector{Float64}`: Masses of the particles.
- `distance_data::Vector{Float64}`: Radial distances of the particles. 
- `z_data::Vector{Float64}`: Metal content of the particles in mass units.
- `max_radius::Float64`: Maximum distance up to which the profile will be calculated.
- `bins::Int64`: Number of subdivisions of [0, `max_radius`] to be used for the profile.

# Returns
- A Tuple of two Arrays.
  The first with the radial distances and the second with the metallicities.
"""
function metallicityProfile(
    mass_data::Vector{Float64},
    distance_data::Vector{Float64},
    z_data::Vector{Float64},
    max_radius::Float64,
    bins::Int64,
)::NTuple{2, Vector{Float64}}

    # Dimension consistency check.
    (
        length(mass_data) == length(distance_data) == length(z_data) ||
        throw(DimensionMismatch("The input vectors should have the same length."))
    )

    # Width of each spherical shell used to calculate the metallicity.
    width = max_radius / bins

    # Initialize output arrays.
    x_data = Vector{Float64}(undef, bins)
    y_data = Vector{Float64}(undef, bins)

    @inbounds for i in eachindex(x_data, y_data)
        # Find the indices of `distance_data` which fall within window `i`.
        idx = findall(x -> width * (i - 1) <= x < width * i, distance_data)

        if isempty(idx)
            x_data[i] = width * (i - 0.5)
            y_data[i] = 0.0
        else
            total_mass = sum(mass_data[idx])
            total_z = sum(z_data[idx])

            # Mean distance for window i.
            x_data[i] = sum(distance_data[idx]) / length(idx)
            # Metallicity for window i.
            y_data[i] = (total_z / total_mass) / SOLAR_METALLICITY
        end
    end

    return x_data, y_data
end

"""
    massProfile(
        mass_data::Vector{Float64},
        distance_data::Vector{Float64},
        max_radius::Float64,
        bins::Int64,
    )::NTuple{2, Vector{Float64}}
	
Compute an accumulated mass profile up to a radius `max_radius`. 

`max_radius` and `distance_data` must be in the same units.

# Arguments
- `mass_data::Vector{Float64}`: Masses of the particles.
- `distance_data::Vector{Float64}`: Radial distances of the particles. 
- `max_radius::Float64`: Maximum distance up to which the profile will be calculated.
- `bins::Int64`: Number of subdivisions of [0, `max_radius`] to be used for the profile.

# Returns
- A Tuple of two Arrays.
  The first with the radial distances and the second with the accumulated masses.
"""
function massProfile(
    mass_data::Vector{Float64},
    distance_data::Vector{Float64},
    max_radius::Float64,
    bins::Int64,
)::NTuple{2, Vector{Float64}}

    # Dimension consistency check.
    (
        length(mass_data) == length(distance_data) ||
        throw(DimensionMismatch("The input vectors should have the same length."))
    )

    # Width of each spherical shell used to calculate the mass.
    width = max_radius / bins

    # Initialize output arrays.
    x_data = Vector{Float64}(undef, bins)
    y_data = Vector{Float64}(undef, bins)

    @inbounds for i in eachindex(x_data, y_data)
        # Indices of `distance_data` within window i.
        idx = findall(x -> width * (i - 1) <= x < width * i, distance_data)
		
		if isempty(idx)
            x_data[i] = width * (i - 0.5)
        else
			# Mean distance for window i.
			x_data[i] = sum(distance_data[idx]) / length(idx)
		end
		
		# Mass for window i.
		y_data[i] = sum(mass_data[idx])
    end

    return x_data, cumsum(y_data)
end

"""
    CMDF(
        mass_data::Vector{Float64},
        metallicity_data::Vector{Float64},
        max_Z::Float64,
        bins::Int64; 
        <keyword arguments>
    )::NTuple{2, Vector{Float64}}
	
Compute the cumulative metallicity distribution function up to a metallicity `max_Z`. 

`mass_data` and `metallicity_data` must be in the same units.

# Arguments
- `mass_data::Vector{Float64}`: Masses of the particles.
- `metallicity_data::Vector{Float64}`: Metallicities of the particles. 
- `max_Z::Float64`: Maximum metallicity up to which the profile will be calculated.
- `bins::Int64`: Number of subdivisions of [0, `max_Z`] to construct the plot.
- `x_norm::Bool = false`: If the x axis will be normalize to its maximum value. 

# Returns
- A Tuple of two Arrays.
  The first with the metallicities and the second with the accumulated masses.
"""
function CMDF(
    mass_data::Vector{Float64},
    metallicity_data::Vector{Float64},
    max_Z::Float64,
    bins::Int64;
    x_norm::Bool = false,
)::NTuple{2, Vector{Float64}}

    # Dimension consistency check.
    (
        length(mass_data) == length(metallicity_data) ||
        throw(DimensionMismatch("The input vectors should have the same length."))
    )
	
	# Dimensionless metallicity.
	Z = metallicity_data ./ mass_data

    # If required, normalize the x axis.
    if x_norm
        Z = Z ./ max_Z
        # Width of the metallicity bins.
        width = 1 / bins
    else
        # Width of the metallicity bins.
	    width = max_Z / bins  
    end
	
	# Total star mass.
	total_m = sum(mass_data)
	
	# Initialize output arrays.
    x_data = Vector{Float64}(undef, bins)
    y_data = Vector{Float64}(undef, bins)

    @inbounds for i in eachindex(x_data, y_data)
		idx = findall(x -> width * (i - 1) <= x < width * i, Z)
		
		if isempty(idx)
            x_data[i] = width * (i - 0.5)
        else
			# Mean metallicity for window i.
			x_data[i] = sum(Z[idx]) / length(idx)
		end
		
		# Mass fraction for window i.
		y_data[i] = sum(mass_data[idx]) / total_m			
		
	end
	
	return x_data, cumsum(y_data)
end

"""
    KennicuttSchmidtLaw(
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
	
Compute mass area density and the SFR area density for the Kennicutt-Schmidt law. 

`temp_filter` and `temperature_data` must be in the same units, and `age_filter` and 
`age_data` must be in the same units too.

# Arguments
- `gas_mass_data::Vector{Float64}`: Masses of the gas particles.
- `gas_distance_data::Vector{Float64}`: 2D distances of the gas particles. 
- `temperature_data::Vector{Float64}`: Temperatures of the gas particles.
- `star_mass_data::Vector{Float64}`: Masses of the stars.
- `star_distance_data::Vector{Float64}`: 2D distances of the stars.
- `age_data::Vector{Float64}`: Ages of the stars.
- `temp_filter::Float64`: Maximum temperature allowed for the gas particles.
- `age_filter::Float64`: Maximum age allowed for the stars.
- `max_r::Float64`: Maximum distance up to which the parameters will be calculated.
- `bins::Int64 = 50`: Number of subdivisions of [0, `max_r`] to be used. 
  It has to be at least 5.

# Returns
- A dictionary with three entries.
  - Key "RHO" => Logarithm of the area mass densities.
  - Key "SFR" => Logarithm of the SFR area densities.
  - Key "LM" => Linear model given by GLM.jl.
"""
function KennicuttSchmidtLaw(
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

    # Bin size check.
    if bins < 5
        error("You have to use at least 10 bins.")
    end

    # Filter out hot gas particles.
    deleteat!(gas_mass_data, temperature_data .> temp_filter)
    deleteat!(gas_distance_data, temperature_data .> temp_filter)
    filter!(x -> x < temp_filter, temperature_data)

    # Filter out old stars.
    deleteat!(star_mass_data, age_data .> age_filter)
    deleteat!(star_distance_data, age_data .> age_filter)
    filter!(x -> x < age_filter, age_data)

    r_width = max_r / bins

    # Initialize output arrays.
    x_data = Vector{Float64}(undef, bins)
    y_data = Vector{Float64}(undef, bins)
    @inbounds for i in eachindex(x_data, y_data)

		idx_gas = findall(x ->  r_width * (i - 1) <= x < r_width * i, gas_distance_data)
        gas_mass = sum(gas_mass_data[idx_gas])
		# Gas area density for window i.
		x_data[i] = gas_mass / (π * r_width * r_width * (2 * i - 1))

        idx_star = findall(x ->  r_width * (i - 1) <= x < r_width * i, star_distance_data)
        sfr = sum(star_mass_data[idx_star]) / age_filter 
        # SFR area density for window i.
        y_data[i] = sfr / (π * r_width * r_width * (2 * i - 1))		
		
	end

    # Filter out zeros.
    deleteat!(y_data, x_data .<= 0.0)
    filter!(x -> x > 0.0, x_data)
    deleteat!(x_data, y_data .<= 0.0)
    filter!(y -> y > 0.0, y_data)

    # Set logarithmic scaling.
    x_data = log10.(x_data)
    y_data = log10.(y_data)

    # If there are too little data to get a fit return nothing
    if length(x_data) < 5
        return nothing
    end

    # Compute linear fit.
    X = [ones(length(x_data)) x_data]
    linear_model = lm(X, y_data)

    return Dict("RHO" => x_data, "SFR" => y_data, "LM" => linear_model)
end

"""
    error_string(mean::Float64, error::Float64)::String

Give the mean and error in a string with the standar formating.

It follows the traditional rules for error printing, the error with only one significant 
digit, unles such digit is a one, in which case, two significant digits are printed. 
The mean will have a presition such as to match the error.

# Arguments 
- `mean::Float64`: Mean value.
- `error::Float64`: Error value. It must be positive.

# Returns
- A Tuple with the formatted mean and error.

# Examples
```julia-repl
julia> error_string(69.42069, 0.038796)
(69.42, 0.04)

julia> error_string(69.42069, 0.018796)
(69.421, 0.019)

julia> error_string(69.42069, 0.0)
(69.42069, 0.0)

julia> error_string(69.42069, 73.4)
(0.0, 70.0)
```
"""
function error_string(mean::Float64, error::Float64)::NTuple{2, Float64}
    error >= 0.0 || error("The error must be a positive number.")

    if error == 0.0
        round_mean = mean
        round_error = error
    else
        sigdigit_pos = abs(log10(abs(error)))
        extra = 0

        if error < 1.0
            if abs(mean) < error
                round_mean = 0.0
            else
                first_dig = trunc(error * 10^(floor(sigdigit_pos) + 1))
                if first_dig == 1.0
                    extra = 1
                end
                digits = ceil(Int64, sigdigit_pos) + extra
                round_mean = round(mean; digits)
            end
        else
            if abs(mean) < error
                round_mean = 0.0
            else
                first_dig = trunc(error / 10^(floor(sigdigit_pos)))
                if first_dig == 1
                    extra = 1
                end

                sigdigits =
                    ceil(Int64, log10(abs(mean))) - ceil(Int64, sigdigit_pos) + 1 + extra
                round_mean = round(mean; sigdigits)
            end
        end

        round_error = round(error, sigdigits = 1 + extra)
    end

    return round_mean, round_error
end