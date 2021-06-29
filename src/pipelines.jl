############################################################################################
# Pipeline functions
############################################################################################

"""
    scatter_grid_pipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64; 
        <keyword arguments>
    )::Nothing

Save the results of the [`scatter_grid_plot`](@ref) function as one image per snapshot, 
and then generate a GIF and video animating the images. 

# Arguments
- `base_name::String`: Base name of the snapshot files, 
  set in the GADGET variable `SnapshotFileBase`.
- `source_path::String`: Path to the directory containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `anim_name::String`: File name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `output_path::String = "scatter_grid"`: Path to the output directory. The images will 
  be stored in `output_path`/images/ and will be named `base_name`\\_XXX`format` where XXX 
  is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. Its unit doesn't have to be the same 
  as `length_unit`.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".ps", ".svg" and ".png". 
"""
function scatter_grid_pipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64; 
    output_path::String = "scatter_grid",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    step::Int64 = 1,
    box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim = get_snapshot_path(base_name, source_path)

    snap_files = @view sim["snap_files"][1:step:end] 
    snap_numbers = @view sim["numbers"][1:step:end] 

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))

    # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Computing the scatter plots... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots                  
    animation = @animate for (number, snapshot) in zip(snap_numbers, snap_files)

        positions = get_position(
            snapshot; 
            sim_cosmo, 
            filter_function, 
            box_size, 
            length_unit,
        )

        figure = scatter_grid_plot(positions)

        savefig(
            figure, 
            joinpath(img_path, base_name * "_" * number * format),
        )

        next!(prog_bar)

    end

    # Make the GIF
    gif(
        animation, 
        joinpath(output_path, anim_name * ".gif"), 
        fps = frame_rate, 
        show_msg = false
    )

    # Make the video
    make_video(img_path, format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    density_map_pipeline(
        base_name::String,
        source_path::String,
        z_quantity::Union{String, Nothing},
        anim_name::String,
        frame_rate::Int64; 
        <keyword arguments>
    )::Nothing

Save the results of the [`density_map_plot`](@ref) function as one image per snapshot, 
and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base name of the snapshot files, 
  set in the GADGET variable `SnapshotFileBase`.
- `source_path::String`: Path to the directory containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `z_quantity::Union{String, Nothing}`: Quantity to be mapped. The options are:
  * `"Z"`: The metallicity (relative to solar metallicity).
  * `"fmol"`: The fraction of molecular gas.
  * `"fatom"`: The fraction of atomic gas.
  * `nothing`: The density itself will be mapped.
- `anim_name::String`: File name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `output_path::String = "density_map"`: Path to the output directory. The images will 
  be stored in `output_path`/images/ and will be named `base_name`\\_XXX`format` where XXX 
  is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `plane::String = "All"`: Indicates which plane will be plotted. 
  * `"XY"` ⟶ x-y plane alone.
  * `"XZ"` ⟶ x-z plane alone.
  * `"YZ"` ⟶ y-z plane alone.
  * `"All"` ⟶ The three planes in a single 1x3 figure.
- `axes::Bool = false`: If true, the axes passing through (0, 0) are drawn. If false, 
  no axes are drawn.
- `box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. Its unit doesn't have to be the same 
  as `length_unit`.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the 
  GR backend can be used, namely ".pdf", ".ps", ".svg" and ".png". 
"""
function density_map_pipeline(
    base_name::String,
    source_path::String,
    z_quantity::Union{String, Nothing},
    anim_name::String,
    frame_rate::Int64;
    output_path::String = "density_map",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    step::Int64 = 1,
    plane::String = "All",
    axes::Bool = false,
    box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim = get_snapshot_path(base_name, source_path)

    snap_files = @view sim["snap_files"][1:step:end]  
    snap_numbers = @view sim["numbers"][1:step:end]  

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))

    # Generate and save the plots                  
    animation = @animate for (number, snapshot) in zip(snap_numbers, snap_files)

        pos = get_position(snapshot; sim_cosmo, filter_function, box_size, length_unit)
        mass = get_mass(snapshot, "gas"; sim_cosmo, filter_function)

        density_unit = mass["unit"] / length_unit^3
    
        ρ = get_density(snapshot; sim_cosmo, filter_function, density_unit)
        hsml = get_hsml(snapshot; sim_cosmo, filter_function, length_unit)

        if z_quantity == "Z"
            metal = get_metallicity(snapshot, "gas"; sim_cosmo, filter_function)
            z = (metal["Z"] ./ mass["mass"]) ./ SOLAR_METALLICITY
        elseif z_quantity == "fmol"
            z = get_fmol(snapshot; sim_cosmo, filter_function)
        elseif z_quantity == "fatom"
            z = get_fatom(snapshot; sim_cosmo, filter_function)
        else
            z = nothing
        end

        figure = density_map_plot(z, pos, mass, ρ, hsml; plane, axes)
 
        savefig(
            figure, 
            joinpath(img_path, base_name * "_" * number * format),
        )

    end

    # Make the GIF
    gif(
        animation, 
        joinpath(output_path, anim_name * ".gif"), 
        fps = frame_rate, 
        show_msg = false
    )

    # Make the video
    make_video(img_path, format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    star_map_pipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64; 
        <keyword arguments>
    )::Nothing

Save the results of the [`star_map_plot`](@ref) function as one image per snapshot, 
and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base name of the snapshot files, 
  set in the GADGET variable `SnapshotFileBase`.
- `source_path::String`: Path to the directory containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `anim_name::String`: File name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `output_path::String = "star_map"`: Path to the output directory. The images will 
  be stored in `output_path`/images/ and will be named `base_name`\\_XXX`format` where XXX 
  is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `plane::String = "All"`: Indicates which plane will be plotted. 
  * `"XY"` ⟶ x-y plane alone.
  * `"XZ"` ⟶ x-z plane alone.
  * `"YZ"` ⟶ y-z plane alone.
  * `"All"` ⟶ The three planes in a single 1x3 figure.
- `box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. Its unit doesn't have to be the same 
  as `length_unit`.
- `box_factor::Float64 = 1.0`: Multiplicative factor for the plotting region. 
  It will scale `positions["box_size"]` if vacuum boundary conditions were used, and
  it will scale `positions["box_size"] / 2` if periodic boundary conditions were used.
- `axes::Bool = false`: If true, the axes passing through (0, 0) are drawn. If false, 
  no axes are drawn.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the 
  GR backend can be used, namely ".pdf", ".ps", ".svg" and ".png". 
"""
function star_map_pipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64;
    output_path::String = "star_map",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    step::Int64 = 1,
    plane::String = "All",
    box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc,
    box_factor::Float64 = 1.0, 
    axes::Bool = false,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,  
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim = get_snapshot_path(base_name, source_path)

    snap_files = @view sim["snap_files"][1:step:end]
    snap_numbers = @view sim["numbers"][1:step:end]

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))

    # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Computing the star map plots... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots                   
    animation = @animate for (number, snapshot) in zip(snap_numbers, snap_files)

        pos = get_position(snapshot; sim_cosmo, filter_function, box_size, length_unit)

        figure = star_map_plot(pos; plane, box_factor, axes)

        savefig(
            figure, 
            joinpath(img_path, base_name * "_" * number * format),
        )

        next!(prog_bar)

    end

    # Make the GIF
    gif(
        animation, 
        joinpath(output_path, anim_name * ".gif"), 
        fps = frame_rate, 
        show_msg = false,
    )

    # Make the video
    make_video(img_path, format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    gas_star_evolution_pipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64; 
        <keyword arguments>
    )::Nothing

Save the results of [`gas_star_evolution_plot`](@ref) function for the last snapshot as one image and 
generate a GIF and a video animating the whole evolution for all snapshots. 
                                
# Arguments
- `base_name::String`: Base name of the snapshot files, 
  set in the GADGET variable `SnapshotFileBase`.
- `source_path::String`: Path to the directory containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `anim_name::String`: File name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `output_path::String = "gas_star_evolution"`: Path to the output directory. 
  The image will be stored in `output_path`/images/ and will be named 
  `base_name`\\_XXX`format` where XXX is the number of the last snapshot. 
  The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots 
  will be plotted.
- `box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. Its unit doesn't have to be the same 
  as `length_unit`.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr`: Unit of mass/time to 
  be used in the output, all available time and mass units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) 
  can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the 
  GR backend can be used, namely ".pdf", ".ps", ".svg" and ".png". 
"""
function gas_star_evolution_pipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64;
    output_path::String = "gas_star_evolution",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    step::Int64 = 1,
    box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr,
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim = get_snapshot_path(base_name, source_path)
    time_series = get_time_evolution(
        sim["snap_files"]; 
        sim_cosmo, 
        filter_function, 
        time_unit, 
        sfr_unit,
    ) 

    snap_files = @view sim["snap_files"][1:step:end]
    snap_numbers = @view sim["numbers"][1:step:end]

    # Create a directory to save the plots, if it doesn't exist
    temp_path = mkpath(joinpath(output_path, "TEMP"))

    # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Computing the gas-star-evolution plots... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots 
    data_iter = enumerate(zip(snap_numbers, snap_files))                 
    animation = @animate for (i, (number, snapshot)) in data_iter

        positions = get_position(
            snapshot; 
            sim_cosmo, 
            filter_function, 
            box_size, 
            length_unit,
        )

        figure = gas_star_evolution_plot(1 + step * (i - 1), time_series, positions)

        savefig(
            figure, 
            joinpath(temp_path, base_name * "_" * number * format),
        )

        next!(prog_bar)

    end

    # Make the GIF
    gif(
        animation, 
        joinpath(output_path, anim_name * ".gif"), 
        fps = frame_rate, 
        show_msg = false,
    )

    # Make the video
    make_video(temp_path, format, output_path, anim_name, frame_rate)

    # Move the last figure out of the temporary directory
    mv(
        joinpath(temp_path, base_name * "_" * snap_numbers[end] * format),
        joinpath(output_path, anim_name * format),
        force = true,
    )

    # Delete the temporary directory and all its contents
    rm(temp_path, recursive = true)

    return nothing
end

"""
    cmdf_pipeline(
        base_name::String,
        source_path::String; 
        <keyword arguments>
    )::Nothing

Save the results of the [`cmdf_plot`](@ref) function as one image per snapshot, if there are 
stars present. 

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `output_path::String = "CMDF"`: Path to the output directory. The images will be stored 
  in `output_path`/images/ and will be named `base_name`\\_XXX`format` where XXX is the 
  number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `x_norm::Bool = false`: If the x axis will be normalized to its maximum value. 
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the time 
  stamps, all available time units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function cmdf_pipeline(
    base_name::String,
    source_path::String;
    output_path::String = "CMDF",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    step::Int64 = 1,
    x_norm::Bool = false,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    format::String = ".png",
)::Nothing

    # Get the simulation data.
    sim = get_snapshot_path(base_name, source_path)
    time_data = get_time_evolution(sim["snap_files"]; sim_cosmo, filter_function, time_unit)

    snap_files = @view sim["snap_files"][1:step:end] 
    snap_numbers = @view sim["numbers"][1:step:end] 
    times = @view time_data["clock_time"][1:step:end]

    # Create a directory to save the plots, if it doesn't exist.
    img_path = mkpath(joinpath(output_path, "images"))

    # Progress bar.
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Computing the CMDF plots... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots.
    data_iter = zip(times, snap_numbers, snap_files)
    for (t, number, snapshot) in data_iter

        header = read_header(snapshot)
        
        if header.nall[5] != 0
            mass_data = get_mass(snapshot, "stars"; sim_cosmo, filter_function)
            z_data = get_metallicity(snapshot, "stars"; sim_cosmo, filter_function)

            figure = cmdf_plot(
                mass_data, 
                z_data,
                t * time_unit; 
                bins = 50, 
                x_norm,
            )

            savefig( 
                figure, 
                joinpath(img_path, base_name * "_" * number * format),
            )

        end

        next!(prog_bar)

    end

    return nothing
end

"""
    cmdf_pipeline(
        base_name::Vector{String},
        source_path::Vector{String},
        labels::Array{String, 2}; 
        <keyword arguments>
    )::Nothing

Save the results of the [`cmdf_plot`](@ref) function for several simulations as one image 
per snapshot, if there are stars present.

# Arguments
- `base_name::Vector{String}`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::Vector{String}`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `labels::Array{String, 2}`: Labels for the different simulations.
- `output_path::String = "CMDF"`: Path to the output directory. The images will be stored 
  in `output_path`/images/ and will be named `base_name`\\_XXX`format` where XXX is the 
  ordinal of the frame. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots 
  will be plotted.
- `x_norm::Bool = false`: If the x axis will be normalized to its maximum value. 
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the time
  stamps, all available time units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function cmdf_pipeline(
    base_name::Vector{String},
    source_path::Vector{String},
    labels::Array{String, 2};
    output_path::String = "CMDF",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    step::Int64 = 1,
    x_norm::Bool = false,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim_data = [get_snapshot_path(base_name[i], path)["snap_files"]
                for (i, path) in enumerate(source_path)]
    # Time stamps (they should be the same for every dataset)
    time_data = get_time_evolution(sim_data[1]; sim_cosmo, filter_function, time_unit)
    times = @view time_data["clock_time"][1:step:end]

    # Length of the shortest simulation
    min_len = minimum(length.(sim_data))
    # Trim longer simulations so all have the same length, 
    # and change their shape to facilitate processing
    trim_matrix = hcat(map(x -> getindex(x, 1:min_len), sim_data)...)
    snap_files = eachrow(trim_matrix[1:step:end, :])

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))
    
    # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Computing the CMDF plots... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots
    data_iter = enumerate(zip(times, snap_files))
    for (i, (t, snapshots)) in data_iter
        
        headers = read_header.(snapshots)
        num_stars = getindex.(getfield.(headers, :nall), 5)

        if all(num_stars .!= 0)
    
            masses = get_mass.(snapshots, "stars"; sim_cosmo, filter_function)
            metallicities = get_metallicity.(snapshots, "stars"; sim_cosmo, filter_function) 

            figure = cmdf_plot(
                masses, 
                metallicities,
                t * time_unit,
                labels;
                bins = 50, 
                x_norm,
            )

            savefig(
                figure, 
                joinpath(img_path, "frame_" * string(step * (i - 1)) * format),
            )

        end

        next!(prog_bar)

    end

    return nothing
end

"""
    birth_histogram_pipeline(
        base_name::String,
        source_path::String; 
        <keyword arguments>
    )::Nothing

Save the results of the [`birth_histogram_plot`](@ref) function as one image per snapshot, 
if there are stars present. 

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `output_path::String = "birth_histogram"`: Path to the output directory. The images 
  will be stored in `output_path`/images/ and will be named `base_name`\\_XXX`format` where 
  XXX is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function birth_histogram_pipeline(
    base_name::String,
    source_path::String;
    output_path::String = "birth_histogram",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    step::Int64 = 1,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim = get_snapshot_path(base_name, source_path)
    time_data = get_time_evolution(sim["snap_files"]; sim_cosmo, filter_function)

    snap_files = @view sim["snap_files"][1:step:end] 
    snap_numbers = @view sim["numbers"][1:step:end] 

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))

    # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Computing the birth histograms... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots
    data_iter = enumerate(zip(snap_numbers, snap_files))
    for (i, (number, snapshot)) in data_iter

        header = read_header(snapshot)

        if header.nall[5] != 0
            nursery = get_birth_place(
                1 + step * (i - 1), 
                sim["snap_files"], 
                time_data["clock_time"],
                time_data["units"]["time"];
                sim_cosmo, 
                filter_function,
                length_unit, 
            )

            figure = birth_histogram_plot(nursery, bins = 50)

            savefig(
                figure, 
                joinpath(img_path, base_name * "_" * number * format),
            )
        end

        next!(prog_bar)

    end

    return nothing
end

@doc raw"""
    evolution_summary_pipeline(
        base_name::String,
        source_path::String,
        fig_name::String; 
        <keyword arguments>
    )::Nothing

Produce up to three figures summarizing the time evolution of the simulation.

The plotted parameters are the number of particles, the total mass, the baryonic 
fractional mass and the star formation rate (SFR), the first three for gas and stars.
If the simulation is Newtonian, only one figure is produced (parameters vs. time),
but if the simulation is cosmological, three figures are produced (parameters vs. time,
parameters vs. scale factor and parameters vs. redshift).

Args:
- `base_name::String`: Base name of the snapshot files, 
  set in the GADGET variable `SnapshotFileBase`.
- `source_path::String`: Path to the directory containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `fig_name::String`: Base name for the figures. The images will be named
  `fig_name`\_vs\_XXX`format` where XXX is 'time', 'redshift' or 'scale_factor'.
- `output_path::String = "evolution_summary"`: Path to the output directory. The images 
  will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `mass_factor::Int64 = 0`: Numerical exponent to scale the mass, e.g. if `mass_factor` = 10 
  the corresponding axis will be scaled by ``10^{10}``.
- `number_factor::Int64 = 0`: Numerical exponent to scale the number of particles, 
  e.g. if `number_factor` = 4 the corresponding axis will be scaled by ``10^4``.
- `mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun`: Unit of mass to be used in the output, 
  all available mass units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr`: Unit of mass/time to 
  be used in the output, all available time and mass units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) 
  can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function evolution_summary_pipeline(
    base_name::String,
    source_path::String,
    fig_name::String;
    output_path::String = "evolution_summary",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    mass_factor::Int64 = 0,
    number_factor::Int64 = 0,
    mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr,
    format::String = ".png",
)::Nothing

    snap_files = get_snapshot_path(base_name, source_path)

    # Create a directory to save the plots, if it doesn't exist
    mkpath(output_path)

    time_series = get_time_evolution(
                    snap_files["snap_files"]; 
                    sim_cosmo, 
                    filter_function,
                    mass_unit, 
                    time_unit, 
                    sfr_unit,
                )

    # Parameters vs. time
    figure_t = time_series_plot(time_series; mass_factor, number_factor)
    
    savefig(
        figure_t, 
        joinpath(output_path, fig_name * "_vs_time" * format),
    )

    if sim_cosmo == 1

        # Parameters vs. scale factor
        figure_a = scale_factor_series_plot(time_series;mass_factor, number_factor)
        
        savefig(
            figure_a, 
            joinpath(output_path, fig_name * "_vs_scale_factor" * format),
        )

        # Parameters vs. redshift
        figure_z = redshift_series_plot(time_series; mass_factor, number_factor)
        
        savefig(
            figure_z, 
            joinpath(output_path, fig_name * "_vs_redshift" * format),
        )

    end

    return nothing
end

@doc raw"""
    compare_simulations_pipeline(
        base_name::Vector{String},
        source_path::Vector{String},
        labels::Array{String, 2},
        fig_name::String,
        x_quantity::String,
        y_quantity::String; 
        <keyword arguments>
    )::Nothing

Make a figure comparing `y_quantity` vs. `x_quantity` for several simulations.

`x_quantity` and `y_quantity` can be any magnitude used in the [`get_time_evolution`](@ref) 
function, namely:

- "scale_factor"                  
- "redshift"                  
- "clock_time" (Physical time)
- "sfr" (SFR) 			              
- "sfr_prob" (SFR probability - Not normalized) 			           
- "gas_number" (Gas particle number) 	            
- "dm_number" (Dark matter particle number)		               
- "star_number" (Star number)         
- "gas_mass" (Total gas mass)               
- "dm_mass" (Total dark matter mass)	            
- "star_mass" (Total star mass)	
- "gas_density" (Total gas density)	               
- "gas_frac" (Gas fraction)		                
- "dm_frac" (Dark matter fraction)		                
- "star_frac" (Star fraction)	                   
- "gas\_bar\_frac" (Baryonic gas fraction)                  
- "star\_bar\_frac" (Baryonic star fraction)

The numeric values of a quantity can also be saved as a text files for the simulations. 
One column per simulation, one row per sanpshot.

# Arguments
- `base_name::Vector{String}`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::Vector{String}`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `labels::Array{String, 2}`: Labels for the different simulations, e.g. [label1 label2 ...].
- `fig_name::String`: Base name for the figure. The file will be named
  `fig_name`_`y_quantity`_vs_`x_quantity` `format`.
- `x_quantity::String`: Physical magnitude for the x axis. 
- `y_quantity::String`: Physical magnitude for the y axis.
- `output_path::String = "compare_simulations"`: Path to the output directory. The images 
  will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `title::String = ""`: Title for the figure. If an empty string is given, no title is 
  printed.
- `x_factor::Int64 = 0`: Numerical exponent to scale the `x_quantity`, e.g. if `x_factor` = 10 
  the corresponding axis will be scaled by ``10^{10}``. The default is no scaling.
- `y_factor::Int64 = 0`: Numerical exponent to scale the `y_quantity`, e.g. if `y_factor` = 10 
  the corresponding axis will be scaled by ``10^{10}``. The default is no scaling.
- `scale::NTuple{2, Symbol} = (:identity, :identity)`: Scaling to be used for the x 
  and y axes. The two options are:
  * `:identity` ⟹ no scaling.
  * `:log10` ⟹ logarithmic scaling.
- `smooth_data::Bool = false`: If true a smoothing window with no weighs is applied to 
  the y data. If false (the default) no transformation occurs.
- `bins::Int64 = 0`: Number of subdivisions for the smoothing of the data, only relevant if
  `smooth_data = true`. 
- `legend_pos::Symbol = :bottomright`: Position of the legend, e.g. `:topleft`.
- `text_quantity::String = ""`: Name of the quantity to be saved in a text file. 
  Any magnitude used in the [`get_time_evolution`](@ref) function can be used. If left empty no 
  text file will be produced.
- `file_name::String = "results"`: Name of the .dat file that will be generated if 
  `text_quantity` is not an empty string.
- `mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun`: Unit of mass to be used in the output, 
  all available mass units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used, 
  e.g. UnitfulAstro.Msun.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used, 
  e.g. UnitfulAstro.Myr.
- `sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr`: Unit of mass/time to 
  be used in the output, all available time and mass units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) 
  can be used, e.g. UnitfulAstro.Msun/UnitfulAstro.yr.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used 
  in the output, all available length units in [UnitfulAstro.jl](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro.jl](https://github.com/JuliaAstro/UnitfulAstro.jl) 
  can be used.
- `density_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.kpc^3`: Unit of 
  density to be used in the output, all available density units in [UnitfulAstro.jl](https://github.com/PainterQubits/Unitful.jl) and 
  [UnitfulAstro.jl](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png".  
"""
function compare_simulations_pipeline(
    base_name::Vector{String},
    source_path::Vector{String},
    labels::Array{String, 2},
    fig_name::String,
    x_quantity::String,
    y_quantity::String;
    output_path::String = "compare_simulations",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    title::String = "",
    x_factor::Int64 = 0,
    y_factor::Int64 = 0,
    scale::NTuple{2, Symbol} = (:identity, :identity),
    smooth_data::Bool = false, 
    bins::Int64 = 50,
    legend_pos::Symbol = :bottomright,
    text_quantity::String = "",
    file_name::String = "results",
    mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    density_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.kpc^3,
    format::String = ".png",
)::Nothing

    time_series = [get_time_evolution(
                        get_snapshot_path(name, path)["snap_files"]; 
                        sim_cosmo,
                        filter_function, 
                        mass_unit, 
                        time_unit, 
                        sfr_unit,
                        length_unit,
                        density_unit,
                    ) for (name, path) in zip(base_name, source_path)]

    # Create a directory to store the figure if it doesn't exist
    mkpath(output_path)

    # Save data in file
    if !isempty(text_quantity)

        file_path = joinpath(output_path, file_name * ".dat")
        # Names of the columns (one per simulation)
        names = join([" \t " * basename(path)  for path in  source_path])
        # Data to be printed out
        # columns = get.(time_series, text_quantity, 0.0)
        columns = [Float64.(col) for col in get.(time_series, text_quantity, 0.0)]
        # Name of the snapshot files (first column)
        snapshots = [
            basename(snap) 
            for snap in get_snapshot_path(base_name[1], source_path[1])["snap_files"]
        ]
        
        # Unit of the quantity
        if text_quantity == "clock_time"
            unit = string(time_unit)
        elseif text_quantity in ["sfr", "sfr_prob"]
            unit = string(sfr_unit)
        elseif text_quantity in ["gas_mass", "dm_mass", "star_mass"]
            unit = string(mass_unit)
        elseif text_quantity == "gas_density"
            unit = string(density_unit)
        else
            unit = "dimensionless"
        end

        open(file_path, "w") do file

            # Fill the datasets with NaN so all have the same length
            max_length = maximum(length.(columns))
            for (i, l_r) in enumerate(length.(columns))
                if l_r < max_length
                    @inbounds for _ in 1:(max_length - l_r)
                        push!(columns[i], NaN)
                    end
                end
            end

            write(file, time_series[1]["labels"][text_quantity] * " [" * unit * "]\n\n")
            write(file, "snapshot" * names * "\n")
            writedlm(file, zip(snapshots, columns...))

        end

    end

    # `y_quantity` vs. `x_quantity` plot
    figure = compare_simulations_plot(
        time_series,
        x_quantity,
        y_quantity,
        labels;
        title,
        x_factor,
        y_factor,
        scale,
        smooth_data,
        bins,
        legend_pos,
    )

    savefig(
        figure, 
        joinpath(output_path, fig_name * "_" * y_quantity * "_vs_" * x_quantity * format),
    )

    return nothing
end

@doc raw"""
    density_histogram_pipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64; 
        <keyword arguments>
    )::Nothing

Save the results of the [`density_histogram_plot`](@ref) function as one image per snapshot, 
and then generate a GIF and a video animating the images. 

Vertical lines with personalized positions and ticks can be added to the plot. By default
none are drawn.

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `anim_name::String`: File name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `output_path::String = "density_histogram"`: Path to the output directory. The images 
  will be stored in `output_path`/images/ and will be named `base_name`\_XXX`format` where 
  XXX is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `flags::Union{Tuple{Vector{<:Real}, Vector{<:AbstractString}}, Nothing} = nothing`: The first 
  vector in the Tuple has the positions of the vetical lines. The second has the 
  corresponding labels. The positions should be in the correct units of density and take 
  into account `factor`.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `factor::Int64 = 0`: Numerical exponent to scale the density, e.g. if `factor` = 10 
  the x axis will be scaled by ``10^{10}``. The default is no scaling.
- `y_scale::Symbol = :identity`: Scaling to be used for the y axis.
  The two options are:
  * `:identity` ⟶ no scaling.
  * `:log10` ⟶ logarithmic scaling.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `density_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.kpc^3`: Unit of 
  density to be used in the output, all available density units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) 
  can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function density_histogram_pipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64;
    output_path::String = "density_histogram",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    flags::Union{Tuple{Vector{<:Real}, Vector{<:AbstractString}}, Nothing} = nothing,
    step::Int64 = 1,
    factor::Int64 = 0,
    y_scale::Symbol = :identity,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    density_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.kpc^3,
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim = get_snapshot_path(base_name, source_path)
    time_data = get_time_evolution(sim["snap_files"]; sim_cosmo, filter_function, time_unit)

    snap_files = @view sim["snap_files"][1:step:end]
    snap_numbers = @view sim["numbers"][1:step:end]
    times = @view time_data["clock_time"][1:step:end]

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))

    # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Computing the density histograms... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots 
    data_iter = zip(times, snap_numbers, snap_files)
    animation = @animate for (t, number, snapshot) in data_iter

        ρ = get_density(snapshot; sim_cosmo, filter_function, density_unit)

        figure = set_vertical_flags(
            flags, 
            density_histogram_plot(ρ, t * time_unit; factor, y_scale),
        )

        savefig(
            figure, 
            joinpath(img_path, base_name * "_" * number * format),
        )

        next!(prog_bar)

    end

    # Make the GIF
    gif(
        animation, 
        joinpath(output_path, anim_name * ".gif"), 
        fps = frame_rate, 
        show_msg = false,
    )

    # Make the video
    make_video(img_path, format, output_path, anim_name, frame_rate)

    return nothing
end

@doc raw"""
    density_profile_pipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64,
        type::String; 
        <keyword arguments>
    )::Nothing

Save the results of the [`density_profile_plot`](@ref) function as one image per snapshot,
and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `anim_name::String`: File name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `type::String`: Particle type.
  * `"gas"` ⟶ Gas particle. 
  * `"dark_matter"` ⟶ Dark matter particle.
  * `"stars"` ⟶ Star particle.
- `output_path::String = "density_profile"`: Path to the output directory. The images 
  will be stored in `output_path`/images/ and will be named `base_name`\_XXX`format` where 
  XXX is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `scale::Symbol = :identity`: Scaling to be used for the y axis.
  The two options are:
  * `:identity` ⟶ no scaling.
  * `:log10` ⟶ logarithmic scaling.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `bins::Int64 = 100`: Number of subdivisions of the region to be used for the profile.
- `factor::Int64 = 0`: Numerical exponent to scale the density, e.g. if `factor` = 10 
  the y axis will be scaled by ``10^{10}``. The default is no scaling.
- `box_factor::Float64 = 1.0`: Multiplicative factor for the plotting region. 
  It will scale `positions["box_size"]` if vacuum boundary conditions were used, and
  it will scale `positions["box_size"] / 2` if periodic boundary conditions were used.
- `box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. Its unit doesn't have to be the same 
  as `length_unit`.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun`: Unit of mass to be used in the output, 
  all available mass units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function density_profile_pipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64,
    type::String;
    output_path::String = "density_profile",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    scale::Symbol = :identity,
    step::Int64 = 1,
    bins::Int64 = 100,
    factor::Int64 = 0,
    box_factor::Float64 = 1.0,
    box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim = get_snapshot_path(base_name, source_path)
    time_data = get_time_evolution(sim["snap_files"]; sim_cosmo, filter_function, time_unit)

    snap_files = @view sim["snap_files"][1:step:end] 
    snap_numbers = @view sim["numbers"][1:step:end] 
    times = @view time_data["clock_time"][1:step:end] 

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))

    # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Computing the density profiles... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots
    data_iter = zip(times, snap_numbers, snap_files)
    animation = @animate for (t, number, snapshot) in data_iter

        positions = get_position(
            snapshot; 
            sim_cosmo, 
            filter_function, 
            box_size, 
            length_unit,
        )
        mass = get_mass(snapshot, type; sim_cosmo, filter_function, mass_unit)

        figure = density_profile_plot(
            positions,
            mass,
            t * time_unit;
            scale,
            bins,
            factor,
            box_factor,
        )

        savefig(
            figure, 
            joinpath(img_path, base_name * "_" * number * format),
        )

        next!(prog_bar)

    end

    # Make the GIF
    gif(
        animation, 
        joinpath(output_path, anim_name * ".gif"), 
        fps = frame_rate, 
        show_msg = false,
    )

    # Make the video
    make_video(img_path, format, output_path, anim_name, frame_rate)

    return nothing
end

@doc raw"""
    density_profile_pipeline(
        base_name::Vector{String},
        source_path::Vector{String},
        anim_name::String,
        frame_rate::Int64,
        type::String,
        labels::Array{String, 2}; 
        <keyword arguments>
    )::Nothing

Save the results of the [`density_profile_plot`](@ref) function for several simulations as one image 
per snapshot, and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::Vector{String}`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::Vector{String}`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `anim_name::String`: File name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `type::String`: Particle type.
  * `"gas"` ⟶ Gas particle. 
  * `"dark_matter"` ⟶ Dark matter particle.
  * `"stars"` ⟶ Star particle.
- `labels::Array{String, 2}`: Labels for the different simulations.
- `output_path::String = "density_profile"`: Path to the output directory. The images will 
  be stored in `output_path`/images/ and will be named frame\_XXX`format` where XXX is the 
  ordinal of the frame. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `scale::Symbol = :identity`: Scaling to be used for the y axis.
  The two options are:
  * `:identity` ⟶ no scaling.
  * `:log10` ⟶ logarithmic scaling.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `bins::Int64 = 100`: Number of subdivisions of the region to be used for the profile.
- `factor::Int64 = 0`: Numerical exponent to scale the density, e.g. if `factor` = 10 
  the y axis will be scaled by ``10^{10}``. The default is no scaling.
- `box_factor::Float64 = 1.0`: Multiplicative factor for the plotting region. 
  It will scale `positions["box_size"]` if vacuum boundary conditions were used, and
  it will scale `positions["box_size"] / 2` if periodic boundary conditions were used.
- `box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. Its unit doesn't have to be the same 
  as `length_unit`.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun`: Unit of mass to be used in the output, 
  all available mass units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function density_profile_pipeline(
    base_name::Vector{String},
    source_path::Vector{String},
    anim_name::String,
    frame_rate::Int64,
    type::String,
    labels::Array{String, 2};
    output_path::String = "density_profile",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    scale::Symbol = :identity,
    step::Int64 = 1,
    bins::Int64 = 100,
    factor::Int64 = 0,
    box_factor::Float64 = 1.0,
    box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    format::String = ".png",
)::Nothing
    
    # Get the simulation data
    sim_data = [get_snapshot_path(base_name[i], path)["snap_files"] 
                for (i, path) in enumerate(source_path)]

    # Time stamps (they should be the same for every dataset)
    time_data = get_time_evolution(sim_data[1]; sim_cosmo, filter_function, time_unit)
    times = @view time_data["clock_time"][1:step:end]

    # Length of the shortest simulation
    min_len = minimum(length.(sim_data))
    # Trim longer simulations so all have the same length, 
    # and change their shape to facilitate processing
    trim_matrix = hcat(map(x -> getindex(x, 1:min_len), sim_data)...)
    snap_files = eachrow(trim_matrix[1:step:end, :])

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))

    # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Computing the density profiles... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots
    data_iter = enumerate(zip(times, snap_files))
    animation = @animate for (i, (t, snapshots)) in data_iter

        positions = get_position.(
            snapshots; 
            sim_cosmo, 
            filter_function, 
            box_size, 
            length_unit,
        )
        masses = get_mass.(snapshots, type; sim_cosmo, filter_function, mass_unit)

        figure = density_profile_plot(
            positions,
            masses,
            t * time_unit,
            labels;
            scale,
            bins,
            factor,
            box_factor,
        )

        savefig(
            figure, 
            joinpath(img_path, "frame_"  * string(step * (i - 1)) * format),
        )

        next!(prog_bar)

    end

    # Make the GIF
    gif(
        animation, 
        joinpath(output_path, anim_name * ".gif"), 
        fps = frame_rate, 
        show_msg = false,
    )

    # Make the video
    make_video(img_path, format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    metallicity_profile_pipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64,
        type::String; 
        <keyword arguments>
    )::Nothing

Save the results of the [`metallicity_profile_plot`](@ref) function as one image per snapshot,
and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `anim_name::String`: File name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `type::String`: Particle type.
  * `"gas"` ⟶ Gas particle. 
  * `"dark_matter"` ⟶ Dark matter particle.
  * `"stars"` ⟶ Star particle.
- `output_path::String = "metallicity_profile"`: Path to the output directory. The images 
  will be stored in `output_path`/images/ and will be named `base_name`\\_XXX`format` where 
  XXX is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `bins::Int64 = 100`: Number of subdivisions of the region to be used for the profile.
- `box_factor::Float64 = 1.0`: Multiplicative factor for the plotting region. 
  It will scale `positions["box_size"]` if vacuum boundary conditions were used, and
  it will scale `positions["box_size"] / 2` if periodic boundary conditions were used.
- `box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. Its unit doesn't have to be the same 
  as `length_unit`.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function metallicity_profile_pipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64,
    type::String;
    output_path::String = "metallicity_profile",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    scale::Symbol = :identity,
    step::Int64 = 1,
    bins::Int64 = 100,
    box_factor::Float64 = 1.0,
    box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim = get_snapshot_path(base_name, source_path)
    time_data = get_time_evolution(sim["snap_files"]; sim_cosmo, filter_function, time_unit)

    snap_files = @view sim["snap_files"][1:step:end] 
    snap_numbers = @view sim["numbers"][1:step:end] 
    times = @view time_data["clock_time"][1:step:end] 

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))

    # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Computing the metallicity profiles... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots
    data_iter = zip(times, snap_numbers, snap_files)
    animation = @animate for (t, number, snapshot) in data_iter

        positions = get_position(
            snapshot; 
            sim_cosmo, 
            filter_function, 
            box_size, 
            length_unit,
        )
        mass = get_mass(snapshot, type; sim_cosmo, filter_function)
        metallicities = get_metallicity(snapshot, type; sim_cosmo, filter_function)

        figure = metallicity_profile_plot(
            positions,
            mass,
            metallicities,
            t * time_unit;
            scale,
            bins,
            box_factor,
        )

        savefig(
            figure, 
            joinpath(img_path, base_name * "_" * number * format),
        )

        next!(prog_bar)

    end

    # Make the GIF
    gif(
        animation, 
        joinpath(output_path, anim_name * ".gif"), 
        fps = frame_rate, 
        show_msg = false,
    )

    # Make the video
    make_video(img_path, format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    metallicity_profile_pipeline(
        base_name::Vector{String},
        source_path::Vector{String},
        anim_name::String,
        frame_rate::Int64,
        type::String,
        labels::Array{String, 2}; 
        <keyword arguments>
    )::Nothing

Save the results of the [`metallicity_profile_plot`](@ref) function for several simulations as one 
image per snapshot, and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::Vector{String}`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::Vector{String}`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `anim_name::String`: File name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `type::String`: Particle type.
  * `"gas"` ⟶ Gas particle. 
  * `"dark_matter"` ⟶ Dark matter particle.
  * `"stars"` ⟶ Star particle.
- `labels::Array{String,2}`: Labels for the different simulations.
- `output_path::String = "metallicity_profile"`: Path to the output directory. The images 
  will be stored in `output_path`/images/ and will be named frame\\_XXX`format` where XXX is 
  the ordinal of the frame. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `scale::Symbol = :identity`: Scaling to be used for the y axis.
  The two options are:
  * `:identity` ⟶ no scaling.
  * `:log10` ⟶ logarithmic scaling.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `bins::Int64 = 100`: Number of subdivisions of the region to be used for the profile.
- `box_factor::Float64 = 1.0`: Multiplicative factor for the plotting region. 
  It will scale `positions["box_size"]` if vacuum boundary conditions were used, and
  it will scale `positions["box_size"] / 2` if periodic boundary conditions were used.
- `box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used.  Its unit doesn't have to be the same 
  as `length_unit`.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function metallicity_profile_pipeline(
    base_name::Vector{String},
    source_path::Vector{String},
    anim_name::String,
    frame_rate::Int64,
    type::String,
    labels::Array{String, 2};
    output_path::String = "metallicity_profile",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    scale::Symbol = :identity,
    step::Int64 = 1,
    bins::Int64 = 100,
    box_factor::Float64 = 1.0,
    box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,  
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim_data = [get_snapshot_path(base_name[i], path)["snap_files"] 
                for (i, path) in enumerate(source_path)]

    # Time stamps (they should be the same for every dataset)
    time_data = get_time_evolution(sim_data[1]; sim_cosmo, filter_function, time_unit)
    times = @view time_data["clock_time"][1:step:end]

    # Length of the shortest simulation
    min_len = minimum(length.(sim_data))
    # Trim longer simulations so all have the same length, 
    # and change their shape to facilitate processing
    trim_matrix = hcat(map(x -> getindex(x, 1:min_len), sim_data)...)
    snap_files = eachrow(trim_matrix[1:step:end, :])

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))

    # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Computing the metallicity profiles... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots
    data_iter = enumerate(zip(times, snap_files))
    animation = @animate for (i, (t, snapshots)) in data_iter

        positions = get_position.(
            snapshots; 
            sim_cosmo, 
            filter_function, 
            box_size, 
            length_unit,
        )
        masses = get_mass.(snapshots, type; sim_cosmo, filter_function)
        metallicities = get_metallicity.(snapshots, type; sim_cosmo, filter_function)

        figure = metallicity_profile_plot(
            positions,
            masses,
            metallicities,
            t * time_unit,
            labels;
            scale,
            bins,
            box_factor,
        )

        savefig(
            figure, 
            joinpath(img_path, "frame_" * string(step * (i - 1)) * format),
        )

        next!(prog_bar)

    end

    # Make the GIF
    gif(
        animation, 
        joinpath(output_path, anim_name * ".gif"), 
        fps = frame_rate, 
        show_msg = false,
    )

    # Make the video
    make_video(img_path, format, output_path, anim_name, frame_rate)

    return nothing
end

@doc raw"""
    mass_profile_pipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64,
        type::String; 
        <keyword arguments>
    )::Nothing

Save the results of the [`mass_profile_plot`](@ref) function as one image per snapshot,
and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `anim_name::String`: File name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `type::String`: Particle type.
  * `"gas"` ⟶ Gas particle. 
  * `"dark_matter"` ⟶ Dark matter particle.
  * `"stars"` ⟶ Star particle.
- `scale::Symbol = :identity`: Scaling to be used for the y axis.
  The two options are:
  * `:identity` ⟶ no scaling.
  * `:log10` ⟶ logarithmic scaling.
- `output_path::String = "mass_profile"`: Path to the output directory. The images will be 
  stored in `output_path`/images/ and will be named `base_name`\_XXX`format` where XXX is 
  the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `bins::Int64 = 100`: Number of subdivisions of the region to be used for the profile.
- `factor::Int64 = 0`: Numerical exponent to scale the mass, e.g. if `factor` = 10 
  the y axis will be scaled by ``10^{10}``. The default is no scaling.
- `box_factor::Float64 = 1.0`: Multiplicative factor for the plotting region. 
  It will scale `positions["box_size"]` if vacuum boundary conditions were used, and
  it will scale `positions["box_size"] / 2` if periodic boundary conditions were used.
- `box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. Its unit doesn't have to be the same 
  as `length_unit`.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun`: Unit of mass to be used in the output, 
  all available mass units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function mass_profile_pipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64,
    type::String;
    output_path::String = "mass_profile",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    scale::Symbol = :identity,
    step::Int64 = 1,
    bins::Int64 = 100,
    factor::Int64 = 0,
    box_factor::Float64 = 1.0,
    box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim = get_snapshot_path(base_name, source_path)
    time_data = get_time_evolution(sim["snap_files"]; sim_cosmo, filter_function, time_unit)

    snap_files = @view sim["snap_files"][1:step:end] 
    snap_numbers = @view sim["numbers"][1:step:end] 
    times = @view time_data["clock_time"][1:step:end] 

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))

    # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Computing the mass profiles... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots
    data_iter = zip(times, snap_numbers, snap_files)
    animation = @animate for (t, number, snapshot) in data_iter

        positions = get_position(
            snapshot; 
            sim_cosmo, 
            filter_function, 
            box_size, 
            length_unit,
        )
        mass = get_mass(snapshot, type; sim_cosmo, filter_function, mass_unit)

        figure = mass_profile_plot(
            positions,
            mass,
            t * time_unit;
            scale,
            bins,
            factor,
            box_factor,
        )

        savefig(
            figure, 
            joinpath(img_path, base_name * "_" * number * format),
        )

        next!(prog_bar)

    end

    # Make the GIF
    gif(
        animation, 
        joinpath(output_path, anim_name * ".gif"), 
        fps = frame_rate, 
        show_msg = false,
    )

    # Make the video
    make_video(img_path, format, output_path, anim_name, frame_rate)

    return nothing
end

@doc raw"""
    mass_profile_pipeline(
        base_name::Vector{String},
        source_path::Vector{String},
        anim_name::String,
        frame_rate::Int64,
        type::String,
        labels::Array{String, 2}; 
        <keyword arguments>
    )::Nothing

Save the results of the [`mass_profile_plot`](@ref) function for several simulations as one image 
per snapshot, and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::Vector{String}`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::Vector{String}`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `anim_name::String`: File name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `type::String`: Particle type.
  * `"gas"` ⟶ Gas particle. 
  * `"dark_matter"` ⟶ Dark matter particle.
  * `"stars"` ⟶ Star particle.
- `labels::Array{String, 2}`: Labels for the different simulations.
- `output_path::String = "mass_profile"`: Path to the output directory. The images will be 
  stored in `output_path`/images/ and will be named frame\_XXX`format` where XXX is the 
  ordinal of the frame. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `scale::Symbol = :identity`: Scaling to be used for the y axis.
  The two options are:
  * `:identity` ⟶ no scaling.
  * `:log10` ⟶ logarithmic scaling.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `bins::Int64 = 100`: Number of subdivisions of the region to be used for the profile.
- `factor::Int64 = 0`: Numerical exponent to scale the mass, e.g. if `factor` = 10 
  the y axis will be scaled by ``10^{10}``. The default is no scaling.
- `box_factor::Float64 = 1.0`: Multiplicative factor for the plotting region. 
  It will scale `positions["box_size"]` if vacuum boundary conditions were used, and
  it will scale `positions["box_size"] / 2` if periodic boundary conditions were used.
- `box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. Its unit doesn't have to be the same 
  as `length_unit`.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun`: Unit of mass to be used in the output, 
  all available mass units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png".  
"""
function mass_profile_pipeline(
    base_name::Vector{String},
    source_path::Vector{String},
    anim_name::String,
    frame_rate::Int64,
    type::String,
    labels::Array{String, 2};
    output_path::String = "mass_profile",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    scale::Symbol = :identity,
    step::Int64 = 1,
    bins::Int64 = 100,
    factor::Int64 = 0,
    box_factor::Float64 = 1.0,
    box_size::Unitful.Quantity = 1000.0UnitfulAstro.kpc,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim_data = [get_snapshot_path(base_name[i], path)["snap_files"] 
                for (i, path) in enumerate(source_path)]

    # Time stamps (they should be the same for every dataset)
    time_data = get_time_evolution(sim_data[1]; sim_cosmo, filter_function, time_unit)
    times = @view time_data["clock_time"][1:step:end]

    # Length of the shortest simulation
    min_len = minimum(length.(sim_data))
    # Trim longer simulations so all have the same length, 
    # and change their shape to facilitate processing
    trim_matrix = hcat(map(x -> getindex(x, 1:min_len), sim_data)...)
    snap_files = eachrow(trim_matrix[1:step:end, :])

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))

    # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Computing the mass profiles... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots
    data_iter = enumerate(zip(times, snap_files))
    animation = @animate for (i, (t, snapshots)) in data_iter

        positions = get_position.(
            snapshots; 
            sim_cosmo, 
            filter_function, 
            box_size, 
            length_unit,
        )
        masses = get_mass.(snapshots, type; sim_cosmo, filter_function)

        figure = mass_profile_plot(
            positions,
            masses,
            t * time_unit,
            labels;
            scale,
            bins,
            factor,
            box_factor,
        )

        savefig(
            figure, 
            joinpath(img_path, "frame_" * string(step * (i - 1)) * format),
        )

        next!(prog_bar)

    end

    # Make the GIF
    gif(
        animation, 
        joinpath(output_path, anim_name * ".gif"), 
        fps = frame_rate, 
        show_msg = false,
    )

    # Make the video
    make_video(img_path, format, output_path, anim_name, frame_rate)

    return nothing
end

@doc raw"""
    sfr_txt_pipeline(
        snapshots::Vector{String},
        source_path::Vector{String},
        x_axis::Int64,
        y_axis::Vector{Int64}; 
        <keyword arguments>
    )::Nothing

Save the results of the [`sfr_txt_plot`](@ref) function as one image per simulation or one image 
per column depending on `comparison_type`.

!!! warning
    This function takes a modified version of sfr.txt which is produced by a private version of 
    GADGET3. GADGET4 produces a sfr.txt, but it is not compatible with this function.

# Arguments
- `snapshots::Vector{String}`: Path to the snapshot files, to get its headers.
- `source_path::Vector{String}`: Paths to the directories containing the sfr.txt files, 
  set in the GADGET variable `OutputDir`.
- `x_axis::Int64`: Column number for the x axis.
- `y_axis::Vector{Int64}`: Column numbers for the y axis.
- `output_path::String = "sfr_txt"`: Path to the output directory.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `comparison_type::Int64 = 0`: Selects which parameters (columns or simulations)
  will be compared:
  * `0` ⟶ Different columns are compared, for a single simulation (one plot per simulation).
  * `1` ⟶ Different simulations are compared, using the same column (one plot per column).
- `titles::Vector{String} = String[]`: Titles for the figures. If an empty string is given,
  no title is printed.
- `names::Vector{String} = String[]`: Names for the files. If an empty string is given, the 
  images will be assigned a number given by the order of `source_path`.
- `labels::Union{Nothing, Array{String, 2}} = nothing`: Labels for the different 
  simulations. Only relevant if `comparison_type = 1`.
- `bins::Int64 = 0`: Number of subdivisions for the smoothing of the data. 
  The default is no smoothing. It will apply equally to every figure produced.
- `scale::NTuple{2, Symbol} = (:identity, :identity)`: Scaling to be used for the x and y 
  axes. It will apply equally to every figure produced.
  The options are:
  * `:identity` ⟶ no scaling.
  * `:log10` ⟶ logarithmic scaling.
- `x_factor::Int64 = 0`: Numerical exponent to scale the `x_quantity`, e.g. if `x_factor` = 10 
  the corresponding axis will be scaled by ``10^{10}``. The default is no scaling.
- `y_factor::Int64 = 0`: Numerical exponent to scale the `y_quantity`, e.g. if `y_factor` = 10 
  the corresponding axis will be scaled by ``10^{10}``. The default is no scaling.
- `min_filter::NTuple{2, Float64} = (-Inf, -Inf)`: Value filter for the x and y axes. 
  It will apply equally to every figure produced. If a value of the x data is lower 
  than `min_filter[1]`, then it is deleted. Equivalently with the y axis and `min_filter[2]`. 
  The default is -Inf for both, i.e. no filtering.
- `mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun`: Unit of mass to be used in the output, 
  all available mass units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr`: Unit of mass/time to 
  be used in the output, all available time and mass units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) 
  can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function sfr_txt_pipeline(
    snapshots::Vector{String},
    source_path::Vector{String},
    x_axis::Int64,
    y_axis::Vector{Int64};
    output_path::String = "sfr_txt",
    sim_cosmo::Int64 = 0,
    comparison_type::Int64 = 0,
    titles::Vector{String} = String[],
    names::Vector{String} = String[],
    labels::Union{Nothing, Array{String, 2}} = nothing,
    bins::Int64 = 0,
    scale::NTuple{2, Symbol} = (:identity, :identity),
    x_factor::Int64 = 0,
    y_factor::Int64 = 0,
    min_filter::NTuple{2, Float64} = (-Inf, -Inf),
    mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr,
    format::String = ".png",
)::Nothing
 
    # By default every figure has no title
    if length(titles) < length(source_path)
        for _ in 1:(length(source_path) - length(titles))
            push!(titles, "")
        end
    end

    # By default every figure has a number as its name
    if length(names) < length(source_path)
        for i in 1:(length(source_path) - length(names))
            push!(names, string(i - 1))
        end
    end

    # Create a directory to save the plots, if it doesn't exist
    mkpath(output_path)

    sfr_data = [
        get_sfr_txt(source, snap; sim_cosmo, mass_unit, time_unit, sfr_unit) 
        for (source, snap) in zip(source_path, snapshots)
    ]

    if comparison_type == 0
        for (data, title, name) in zip(sfr_data, titles, names)

            figure = sfr_txt_plot(
                data, 
                x_axis, 
                y_axis;
                title, 
                bins, 
                scale,
                x_factor, 
                y_factor, 
                min_filter,
            )
            
            savefig(
                figure, 
                joinpath(output_path, name * format),
            )

        end
    else 
        @inbounds for (column, title, name) in zip(y_axis, titles, names)

            figure = sfr_txt_plot(
                sfr_data, 
                x_axis, 
                column,
                labels;
                title, 
                bins, 
                scale,
                x_factor, 
                y_factor, 
                min_filter,
            )

            savefig(
                figure, 
                joinpath(output_path, name * format),
            )

        end
    end

    return nothing
end

"""
    temperature_histogram_pipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64; 
        <keyword arguments>
    )::Nothing

Save the results of the [`temperature_histogram_plot`](@ref) function as one image per snapshot, 
and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `anim_name::String`: File name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `output_path::String = "temperature_histogram"`: Path to the output directory. The images 
  will be stored in `output_path`/images/ and will be named `base_name`\\_XXX`format` where 
  XXX is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `temp_unit::Unitful.FreeUnits = Unitful.K`: Unit of temperature to be used in the 
  output, all available temperature units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function temperature_histogram_pipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64;
    output_path::String = "temperature_histogram",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    step::Int64 = 1,
    temp_unit::Unitful.FreeUnits = Unitful.K,
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim = get_snapshot_path(base_name, source_path)
    time_data = get_time_evolution(sim["snap_files"]; sim_cosmo, filter_function)
    time_unit = time_data["units"]["time"]

    snap_files = @view sim["snap_files"][1:step:end] 
    snap_numbers = @view sim["numbers"][1:step:end] 
    times = @view time_data["clock_time"][1:step:end]

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))
    
    # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Computing the temperature histograms... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots
    data_iter = zip(times, snap_numbers, snap_files)
    animation = @animate for (t, number, snapshot) in data_iter

        temp_data = get_temperature(snapshot; sim_cosmo, filter_function, temp_unit)

        figure = temperature_histogram_plot(temp_data, t * time_unit, bins = 30)

        savefig(
            figure, 
            joinpath(img_path, base_name * "_" * number * format),
        )

        next!(prog_bar)
        
    end

    # Make the GIF
    gif(
        animation, 
        joinpath(output_path, anim_name * ".gif"), 
        fps = frame_rate, 
        show_msg = false,
    )

    # Make the video
    make_video(img_path, format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    rho_temp_pipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64; 
        <keyword arguments>
    )::Nothing

Save the results of the [`rho_temp_plot`](@ref) function as one image per snapshot, 
and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `anim_name::String`: File name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `output_path::String = "rho_vs_temp"`: Path to the output directory. The images will be 
  stored in `output_path`/images/ and will be named `base_name`\\_XXX`format` where XXX is 
  the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `temp_unit::Unitful.FreeUnits = Unitful.K`: Unit of temperature to be used in the 
  output, all available temperature units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `density_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.kpc^3`: Unit of 
  density to be used in the output, all available density units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) 
  can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function rho_temp_pipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64;
    output_path::String = "rho_vs_temp",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    step::Int64 = 1,
    temp_unit::Unitful.FreeUnits = Unitful.K,
    density_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.kpc^3,
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim = get_snapshot_path(base_name, source_path)
    time_data = get_time_evolution(sim["snap_files"]; sim_cosmo, filter_function)
    time_unit = time_data["units"]["time"]

    snap_files = @view sim["snap_files"][1:step:end] 
    snap_numbers = @view sim["numbers"][1:step:end] 
    times = @view time_data["clock_time"][1:step:end]

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))
    
    # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Generating the ρ vs T plots... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots
    data_iter = zip(times, snap_numbers, snap_files)
    animation = @animate for (t, number, snapshot) in data_iter

        temp_data = get_temperature(snapshot; sim_cosmo, filter_function, temp_unit)
        density_data = get_density(snapshot; sim_cosmo, filter_function, density_unit)

        figure = rho_temp_plot(temp_data, density_data, t * time_unit)

        savefig(
            figure, 
            joinpath(img_path, base_name * "_" * number * format),
        )

        next!(prog_bar)
        
    end

    # Make the GIF
    gif(
        animation, 
        joinpath(output_path, anim_name * ".gif"), 
        fps = frame_rate, 
        show_msg = false,
    )

    # Make the video
    make_video(img_path, format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    fraction_temp_pipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64,
        fraction::String; 
        <keyword arguments>
    )::Nothing

Save the results of the [`fraction_temp_plot`](@ref) function as one image per snapshot, 
and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `anim_name::String`: File name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `fraction::String`: Which fraction to plot against the temperature. The options are
  * `"atomic"`: Atomic fraction.
  * `"molecular"`: Molecular fraction.
- `output_path::String = "rho_vs_temp"`: Path to the output directory. The images will be 
  stored in `output_path`/images/ and will be named `base_name`\\_XXX`format` where XXX is 
  the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `temp_unit::Unitful.FreeUnits = Unitful.K`: Unit of temperature to be used in the 
  output, all available temperature units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function fraction_temp_pipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64,
    fraction::String;
    output_path::String = "rho_vs_temp",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    step::Int64 = 1,
    temp_unit::Unitful.FreeUnits = Unitful.K,
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim = get_snapshot_path(base_name, source_path)
    time_data = get_time_evolution(sim["snap_files"]; sim_cosmo, filter_function)
    time_unit = time_data["units"]["time"]

    snap_files = @view sim["snap_files"][1:step:end] 
    snap_numbers = @view sim["numbers"][1:step:end] 
    times = @view time_data["clock_time"][1:step:end]

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))
    
    # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Generating the fmol/fatom vs T plots... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots
    data_iter = zip(times, snap_numbers, snap_files)
    animation = @animate for (t, number, snapshot) in data_iter

        temp_data = get_temperature(snapshot; sim_cosmo, filter_function, temp_unit)

        if fraction == "molecular"
            fmol = get_fmol(snapshot; sim_cosmo, filter_function)
            figure = fraction_temp_plot(temp_data, fmol, t * time_unit, "Molecular fraction")
        elseif fraction == "atomic"
            fatom = get_fmol(snapshot; sim_cosmo, filter_function)
            figure = fraction_temp_plot(temp_data, fatom, t * time_unit, "Atomic fraction")
        end

        savefig(
            figure, 
            joinpath(img_path, base_name * "_" * number * format),
        )

        next!(prog_bar)
        
    end

    # Make the GIF
    gif(
        animation, 
        joinpath(output_path, anim_name * ".gif"), 
        fps = frame_rate, 
        show_msg = false,
    )

    # Make the video
    make_video(img_path, format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    kennicutt_schmidt_pipeline(
        base_name::String,
        source_path::String; 
        <keyword arguments>
    )::Nothing

Save the results of the [`kennicutt_schmidt_plot`](@ref) function as one image per snapshot.

It will produce output only for the snapshots that have enough young stars to produce 
at least five data points for the linear fitting.

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `output_path::String = "Kennicutt_Schmidt"`: Path to the output directory. The images will 
  be stored in `output_path`/images/ and will be named `base_name`\\_XXX`format` where XXX 
  is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `temp_filter::Unitful.Quantity = Inf * Unitful.K`: Maximum temperature allowed for the 
  gas particles.
- `age_filter::Unitful.Quantity = 20.0UnitfulAstro.Myr`: Maximum star age allowed for the 
  calculation of the SFR.
- `max_r::Unitful.Quantity = 1000.0UnitfulAstro.kpc`: Maximum distance up to which the 
  parameters will be calculated, with units.
- `bins::Int64 = 50`: Number of subdivisions of [0, `max_r`] to be used. 
  It has to be at least 5.
- `error_formating::String = "std_error"`: What to print as error values. The options are:
  * `"std_error"` ⟹ mean ± standard_error.
  * `"conf_interval"` ⟹ mean ± max(upper\\_95% - mean, mean - lower\\_95%).
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun`: Unit of mass to be used in the output, 
  all available mass units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function kennicutt_schmidt_pipeline(
    base_name::String,
    source_path::String;
    output_path::String = "Kennicutt_Schmidt",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    step::Int64 = 1,
    temp_filter::Unitful.Quantity = Inf * Unitful.K,
    age_filter::Unitful.Quantity = 20.0UnitfulAstro.Myr,
    max_r::Unitful.Quantity = 1000.0UnitfulAstro.kpc,
    bins::Int64 = 50,
    error_formating::String = "std_error",
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim = get_snapshot_path(base_name, source_path)
    time_data = get_time_evolution(sim["snap_files"]; sim_cosmo, filter_function, time_unit)

    snap_files = @view sim["snap_files"][1:step:end] 
    snap_numbers = @view sim["numbers"][1:step:end] 
    times = @view time_data["clock_time"][1:step:end]

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))

     # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Generating the Kennicutt-Schmidt plots... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots
    data_iter = zip(times, snap_numbers, snap_files)   
    # Initial snapshot
    snap_0 = first(snap_files)     
    for (t, number, snapshot) in data_iter

        header = read_header(snapshot)
        if header.nall[5] != 0

            # Gas masses
            gas_mass_data = get_mass(snapshot, "gas"; sim_cosmo, filter_function, mass_unit)
            # Gas temperatures
            temperature_data = get_temperature(
                snapshot; 
                sim_cosmo, 
                filter_function, 
                temp_unit = unit(temp_filter),
            )
            # Stars masses
            star_mass_data = get_mass(
                snapshot, 
                "stars"; 
                sim_cosmo, 
                filter_function, 
                mass_unit,
            )
            # Stars ages
            age_data = get_age(
                snapshot, 
                t * time_unit; 
                sim_cosmo, 
                snap_0, 
                filter_function,
            )
            # Positions
            pos_data = get_position(
                snapshot; 
                sim_cosmo, 
                filter_function,  
                length_unit,
            )

            figure = kennicutt_schmidt_plot(
                gas_mass_data,
                temperature_data,
                star_mass_data,
                age_data,
                pos_data,
                temp_filter,
                age_filter,
                max_r,
                t * time_unit;
                bins,
                error_formating,
            )

            if figure !== nothing
                # If there was enough data to make a fit
                
                savefig(
                    figure, 
                    joinpath(img_path, base_name * "_" * number * format),
                )
            end

        end

        next!(prog_bar)

    end

    return nothing
end

"""
    cpu_txt_pipeline(
        source_path::Vector{String},
        targets::Vector{String}; 
        <keyword arguments>
    )::Nothing

Save the result of the [`cpu_txt_plot`](@ref) function as one image per simulation.

# Arguments
- `source_path::Vector{String}`: Paths to the directories containing the cpu.txt files, 
  set in the GADGET variable `OutputDir`.
- `targets::Vector{String}`: Target processes to be plotted for each simulation.
- `step::Int64 = 1`: Step used to traverse the CPU cycles, i.e. one every `step` cycles is 
  used for the output plot.
- `output_path::String = "cpu_txt"`: Path to the output directory.
- `title::Vector{String} = String[]`: Titles for the figures. If an empty string is given 
  no title is printed, which is the default.
- `names::Vector{String} = String[]`: Names for the files. If an empty string is given, the 
  images will be assigned a number, starting from 0.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function cpu_txt_pipeline(
    source_path::Vector{String},
    targets::Vector{String};
    step::Int64 = 1,
    output_path::String = "cpu_txt",
    titles::Vector{String} = String[],
    names::Vector{String} = String[],
    format::String = ".png",
)::Nothing

    # By default every figure has no title
    if length(titles) < length(source_path)
        for _ in 1:(length(source_path) - length(titles))
            push!(titles, "")
        end
    end

    # By default every figure has a number as its name
    if length(names) < length(source_path)
        for i in 1:(length(source_path) - length(names))
            push!(names, string(i - 1))
        end
    end

    # Create a directory to save the plots, if it doesn't exist
    mkpath(output_path)

    # Generate and save the plots
    for (source, title, name) in zip(source_path, titles, names)

        figure = cpu_txt_plot(get_cpu_txt(source, targets; step) , title)

        
        savefig(
            figure, 
            joinpath(output_path, name * format),
        )

    end

    return nothing
end

"""
    cpu_txt_pipeline(
        source_path::Vector{String},
        target::String,
        labels::Array{String, 2}; 
        <keyword arguments>
    )::Nothing

Save the result of the [`cpu_txt_plot`](@ref) function, comparing the CPU usage of one 
process among several simulations. 

# Arguments
- `source_path::Vector{String}`: Paths to the directories containing the cpu.txt files, 
  set in the GADGET variable `OutputDir`.
- `target::String`: Target process.
- `labels::Array{String, 2}`: Labels for the different simulations.
- `step::Int64 = 1`: Step used to traverse the CPU cycles, i.e. one every `step` cycles is 
  used for the output plot.
- `output_path::String = "cpu_txt"`: Path to the output directory.
- `title::String = ""`: Title for the figure. If an empty string is given no title is 
  printed, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function cpu_txt_pipeline(
    source_path::Vector{String},
    target::String,
    labels::Array{String, 2};
    step::Int64 = 1,
    output_path::String = "cpu_txt",
    title::String = "",
    format::String = ".png",
)::Nothing

    # Create a directory to save the plots, if it doesn't exist
    mkpath(output_path)
    
    data = [get_cpu_txt(source, target; step) for source in source_path]

    # Generate and save the plots
    figure = cpu_txt_plot(data, labels, title)

    savefig(
        figure, 
        joinpath(output_path, "compare_cpu_txt" * format),
    )

    return nothing
end

"""
    quantities_2D_pipeline(
        base_name::String,
        source_path::String; 
        <keyword arguments>
    )

Save the results of the [`quantities_2D_plot`](@ref) function as one folder per snapshot.

It will produce output only for snapshots that have stars.

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `output_path::String = "Kennicutt_Schmidt"`: Path to the output directory. The images will 
  be stored in `output_path`/images/ and will be named `base_name`\\_XXX`format` where XXX 
  is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `title::String = ""`: Title for the figure. If an empty string is given no title is 
  printed, which is the default.
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `temp_filter::Unitful.Quantity = Inf * Unitful.K`: Maximum temperature allowed for the 
  gas particles.
- `age_filter::Unitful.Quantity = 20.0UnitfulAstro.Myr`: Maximum star age allowed for the 
  calculation of the SFR.
- `max_r::Unitful.Quantity = 1000.0UnitfulAstro.kpc`: Maximum distance up to which the 
  parameters will be calculated, with units.
- `bins::Int64 = 50`: Number of subdivisions of [0, `max_r`] to be used. 
  It has to be at least 5.
- `scale::NTuple{2, Symbol} = (:identity, :identity)`: Scaling to be used for the x and y 
  axes. It will apply equally to every figure produced.
  The options are:
  * `:identity` ⟶ no scaling.
  * `:log10` ⟶ logarithmic scaling.
- `x_factor::Int64 = 0`: Numerical exponent to scale the `x_quantity`, e.g. if `x_factor` = 10 
  the corresponding axis will be scaled by ``10^{10}``. The default is no scaling.
- `y_factor::Int64 = 0`: Numerical exponent to scale the `y_quantity`, e.g. if `y_factor` = 10 
  the corresponding axis will be scaled by ``10^{10}``. The default is no scaling.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun`: Unit of mass to be used in the output, 
  all available mass units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function quantities_2D_pipeline(
    base_name::String,
    source_path::String;
    output_path::String = "quantities_2D",
    sim_cosmo::Int64 = 0,
    title::String = "",
    filter_function::Function = pass_all,
    step::Int64 = 1,
    temp_filter::Unitful.Quantity = Inf * Unitful.K,
    age_filter::Unitful.Quantity = 20.0UnitfulAstro.Myr,
    max_r::Unitful.Quantity = 1000.0UnitfulAstro.kpc,
    bins::Int64 = 50,
    scale::NTuple{2, Symbol} = (:identity, :identity),
    x_factor::Int64 = 0,
    y_factor::Int64 = 0,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    format::String = ".png",
)::Nothing

    x_quantities = ["STARS", "P"]
    y_quantities = ["SFR", "SSFR", "SFE", "GAS", "Psi_FMOL", "OH"]

    # Get the simulation data
    sim = get_snapshot_path(base_name, source_path)
    time_data = get_time_evolution(sim["snap_files"]; sim_cosmo, filter_function, time_unit)

    snap_files = @view sim["snap_files"][1:step:end] 
    snap_numbers = @view sim["numbers"][1:step:end] 
    times = @view time_data["clock_time"][1:step:end]

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(output_path)

     # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Generating the quantities-2D plots... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots
    data_iter = zip(times, snap_numbers, snap_files)   
    # Initial snapshot
    snap_0 = first(snap_files)     
    for (t, number, snapshot) in data_iter

        header = read_header(snapshot)
        if header.nall[5] != 0

            # Gas masses
            gas_mass_data = get_mass(snapshot, "gas"; sim_cosmo, filter_function, mass_unit)
            # Gas temperatures
            temperature_data = get_temperature(
                snapshot; 
                sim_cosmo, 
                filter_function, 
                temp_unit = unit(temp_filter),
            )
            # Stars masses
            star_mass_data = get_mass(
                snapshot, 
                "stars"; 
                sim_cosmo, 
                filter_function, 
                mass_unit,
            )
            # Stars ages
            age_data = get_age(
                snapshot, 
                t * time_unit; 
                sim_cosmo, 
                snap_0, 
                filter_function,
            )
            # Distances
            pos_data = get_position(
                snapshot; 
                sim_cosmo, 
                filter_function,  
                length_unit,
            )
            gas_distance = [norm(col) for col in eachcol(pos_data["gas"])]
            star_distance =  [norm(col) for col in eachcol(pos_data["stars"])]
            # Molecular fraction
            fmol = get_fmol(
                snapshot; 
                sim_cosmo, 
                filter_function,
            )
            # Gas metal mass
            gas_mz = get_metal_mass(
                snapshot, 
                "gas";
                sim_cosmo,
                filter_function,
                mass_unit,
            )

            quantities2D = GADGETPlotting.quantities_2D(
                gas_mass_data["mass"],
                gas_distance,
                temperature_data["temperature"],
                star_mass_data["mass"],
                star_distance,
                age_data["ages"],
                gas_mz["Z"],
                fmol,
                ustrip(temp_filter),
                ustrip(age_filter),	
                ustrip(max_r);
                bins,
            )

            snap_folder = mkpath(joinpath(img_path, base_name * "_" * number))

            for x_quantitie in x_quantities
                for y_quantitie in y_quantities

                    figure = quantities_2D_plot(
                        quantities2D,
                        x_quantitie,
                        y_quantitie,
                        Dict(
                            "mass" => mass_unit, 
                            "length" => length_unit, 
                            "time" => time_unit,
                        );
                        title,
                        x_factor,
                        y_factor,
                        scale,
                    )
                    savefig(
                        figure, 
                        joinpath(
                            snap_folder, 
                            y_quantitie * "_vs_" * x_quantitie * format,
                        ),
                    )

                end
            end

        end

        next!(prog_bar)

    end

    return nothing
end

"""
    fraction_histogram_pipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64,
        fraction::String; 
        <keyword arguments>
    )::Nothing

Save the results of the [`fraction_histogram_plot`](@ref) function as one image per snapshot, 
and then generate a GIF and a video animating the images. 

Vertical lines with personalized positions and ticks can be added to the plot. By default
none are drawn.

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `anim_name::String`: File name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `fraction::String`: Which fraction to plot against the temperature. The options are
  * `"atomic"`: Atomic fraction.
  * `"molecular"`: Molecular fraction.
- `output_path::String = "density_histogram"`: Path to the output directory. The images 
  will be stored in `output_path`/images/ and will be named `base_name`\\_XXX`format` where 
  XXX is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `flags::Union{Tuple{Vector{<:Real}, Vector{<:AbstractString}}, Nothing} = nothing`: The first 
  vector in the Tuple has the positions of the vetical lines. The second has the 
  corresponding labels. The positions should be in the correct units of density and take 
  into account `factor`.
- `bins::Int64 = 20`: Number of subdivisions used for the histogram.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `y_scale::Symbol = :identity`: Scaling to be used for the y axis.
  The two options are:
  * `:identity` ⟶ no scaling.
  * `:log10` ⟶ logarithmic scaling.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in [Unitful](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl) can be used.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function fraction_histogram_pipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64,
    fraction::String;
    output_path::String = "fraction_histogram",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    flags::Union{Tuple{Vector{<:Real}, Vector{<:AbstractString}}, Nothing} = nothing,
    bins::Int64 = 20,
    step::Int64 = 1,
    y_scale::Symbol = :identity,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim = get_snapshot_path(base_name, source_path)
    time_data = get_time_evolution(sim["snap_files"]; sim_cosmo, filter_function, time_unit)

    snap_files = @view sim["snap_files"][1:step:end]
    snap_numbers = @view sim["numbers"][1:step:end]
    times = @view time_data["clock_time"][1:step:end]

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))

    # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Computing the fmol/fatom histograms... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots 
    data_iter = zip(times, snap_numbers, snap_files)
    animation = @animate for (t, number, snapshot) in data_iter

        if fraction == "molecular"

            fmol = get_fmol(snapshot; sim_cosmo, filter_function)
            figure = fraction_histogram_plot(
                fmol, 
                t * time_unit, 
                "Molecular fraction";
                bins, 
                y_scale,
            )

        elseif fraction == "atomic"

            fatom = get_fmol(snapshot; sim_cosmo, filter_function)
            figure = fraction_histogram_plot(
                fatom, 
                t * time_unit, 
                "Atomic fraction"; 
                bins, 
                y_scale,
            )

        end

        figure = set_vertical_flags(flags, figure)

        savefig(
            figure, 
            joinpath(img_path, base_name * "_" * number * format),
        )

        next!(prog_bar)

    end

    # Make the GIF
    gif(
        animation, 
        joinpath(output_path, anim_name * ".gif"), 
        fps = frame_rate, 
        show_msg = false,
    )

    # Make the video
    make_video(img_path, format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    fmol_fatom_pipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64; 
        <keyword arguments>
    )::Nothing

Save the results of the [`fmol_fatom_plot`](@ref) function as one image per snapshot, 
and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `anim_name::String`: File name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `output_path::String = "fmol_vs_ftom"`: Path to the output directory. The images 
  will be stored in `output_path`/images/ and will be named `base_name`\\_XXX`format` where 
  XXX is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png". 
"""
function fmol_fatom_pipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64;
    output_path::String = "fmol_vs_ftom",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    step::Int64 = 1,
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim = get_snapshot_path(base_name, source_path)
    time_data = get_time_evolution(sim["snap_files"]; sim_cosmo, filter_function)
    time_unit = time_data["units"]["time"]

    snap_files = @view sim["snap_files"][1:step:end] 
    snap_numbers = @view sim["numbers"][1:step:end] 
    times = @view time_data["clock_time"][1:step:end]

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))
    
    # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Generating the fmol vs fatom plots... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots
    data_iter = zip(times, snap_numbers, snap_files)
    animation = @animate for (t, number, snapshot) in data_iter

        density = get_density(snapshot; sim_cosmo, filter_function)
        fmol = get_fmol(snapshot; sim_cosmo, filter_function)
        fatom = get_fmol(snapshot; sim_cosmo, filter_function)
        mass = get_mass(snapshot, "gas"; sim_cosmo, filter_function)
        metal = get_metallicity(snapshot, "gas"; sim_cosmo, filter_function)
        
        figure = fmol_fatom_plot(fmol, fatom, density, metal, mass, t * time_unit)

        savefig(
            figure, 
            joinpath(img_path, base_name * "_" * number * format),
        )

        next!(prog_bar)
        
    end

    # Make the GIF
    gif(
        animation, 
        joinpath(output_path, anim_name * ".gif"), 
        fps = frame_rate, 
        show_msg = false,
    )

    # Make the video
    make_video(img_path, format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    fatom_rho_pipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64; 
        <keyword arguments>
    )::Nothing

Save the results of the [`fatom_rho_plot`](@ref) function as one image per snapshot, 
and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `anim_name::String`: File name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `output_path::String = "fatom_vs_rho""`: Path to the output directory. The images 
  will be stored in `output_path`/images/ and will be named `base_name`\\_XXX`format` where 
  XXX is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png".  
"""
function fatom_rho_pipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64;
    output_path::String = "fatom_vs_rho",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    step::Int64 = 1,
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim = get_snapshot_path(base_name, source_path)
    time_data = get_time_evolution(sim["snap_files"]; sim_cosmo, filter_function)
    time_unit = time_data["units"]["time"]

    snap_files = @view sim["snap_files"][1:step:end] 
    snap_numbers = @view sim["numbers"][1:step:end] 
    times = @view time_data["clock_time"][1:step:end]

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))
    
    # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Generating the fatom vs ρ plots... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots
    data_iter = zip(times, snap_numbers, snap_files)
    animation = @animate for (t, number, snapshot) in data_iter

        density = get_density(snapshot; sim_cosmo, filter_function)
        fatom = get_fatom(snapshot; sim_cosmo, filter_function)
        
        figure = fatom_rho_plot(fatom, density, t * time_unit)

        savefig(
            figure, 
            joinpath(img_path, base_name * "_" * number * format),
        )

        next!(prog_bar)
        
    end

    # Make the GIF
    gif(
        animation, 
        joinpath(output_path, anim_name * ".gif"), 
        fps = frame_rate, 
        show_msg = false,
    )

    # Make the video
    make_video(img_path, format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    fmol_Z_pipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64; 
        <keyword arguments>
    )::Nothing

Save the results of the [`fatom_rho_plot`](@ref) function as one image per snapshot, 
and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable `SnapshotFileBase`.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable `OutputDir`.
- `anim_name::String`: File name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `output_path::String = "fmol_vs_Z""`: Path to the output directory. The images 
  will be stored in `output_path`/images/ and will be named `base_name`\\_XXX`format` where 
  XXX is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable `ComovingIntegrationOn`: 
  * `0` ⟶ Newtonian simulation (static universe).
  * `1` ⟶ Cosmological simulation (expanding universe).
- `filter_function::Function = pass_all`: A function with the signature: 

  `foo(snap_file::String, type::String)::Vector{Int64}`
  
  See the function [`pass_all`](@ref) for an example. By default, no particles are filtered.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. By default all snapshots will be plotted.
- `format::String = ".png"`: File format of the output figure. All formats supported by the
  GR backend can be used, namely ".pdf", ".svg" and ".png".  
"""
function fmol_Z_pipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64;
    output_path::String = "fmol_vs_Z",
    sim_cosmo::Int64 = 0,
    filter_function::Function = pass_all,
    step::Int64 = 1,
    format::String = ".png",
)::Nothing

    # Get the simulation data
    sim = get_snapshot_path(base_name, source_path)
    time_data = get_time_evolution(sim["snap_files"]; sim_cosmo, filter_function)
    time_unit = time_data["units"]["time"]

    snap_files = @view sim["snap_files"][1:step:end] 
    snap_numbers = @view sim["numbers"][1:step:end] 
    times = @view time_data["clock_time"][1:step:end]

    # Create a directory to save the plots, if it doesn't exist
    img_path = mkpath(joinpath(output_path, "images"))
    
    # Progress bar
    prog_bar = Progress(
        length(snap_files), 
        dt = 0.5, 
        desc = "Generating the fmol vs Z plots... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
    )

    # Generate and save the plots
    data_iter = zip(times, snap_numbers, snap_files)
    animation = @animate for (t, number, snapshot) in data_iter

        mass = get_mass(snapshot, "gas"; sim_cosmo, filter_function)
        metal = get_metallicity(snapshot, "gas"; sim_cosmo, filter_function)
        fmol = get_fmol(snapshot; sim_cosmo, filter_function)
        
        figure = fmol_Z_plot(fmol, metal, mass, t * time_unit)

        savefig(
            figure, 
            joinpath(img_path, base_name * "_" * number * format),
        )

        next!(prog_bar)
        
    end

    # Make the GIF
    gif(
        animation, 
        joinpath(output_path, anim_name * ".gif"), 
        fps = frame_rate, 
        show_msg = false,
    )

    # Make the video
    make_video(img_path, format, output_path, anim_name, frame_rate)

    return nothing
end