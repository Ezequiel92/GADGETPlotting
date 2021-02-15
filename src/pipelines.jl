############################################################################################
# PIPELINE FUNCTIONS.
############################################################################################

"""
    scatterGridPipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64; 
        <keyword arguments>
    )::Nothing

Save the results of the scatterGridPlot function as one image per snapshot, 
and then generate a GIF and video animating the images. 

# Arguments
- `base_name::String`: Base name of the snapshot files, 
  set in the GADGET variable SnapshotFileBase.
- `source_path::String`: Path to the directory containing the snapshot files, 
  set in the GADGET variable OutputDir.
- `anim_name::String`: Name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `output_path::String = "scatter_grid/"`: Path to the output directory. The images will 
  be stored in `output_path`images/ and will be named `base_name`_XXX`format` where XXX 
  is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `step::Int64=1`: Step used to traverse the list of snapshots. The default is 1, 
  i.e. all snapshots will be plotted.
- `box_size::Unitful.Quantity = 1000UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. It has to have units, e.g. 1000UnitfulAstro.kpc, 
  which is the default. Its units don't have to be the same as `length_unit`.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.kpc, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by GR 
  can be used, namely ".pdf", ".ps", ".svg" and ".png", which is the default. 
"""
function scatterGridPipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64; 
    output_path::String = "scatter_grid/",
    sim_cosmo::Int64 = 0,
    step::Int64 = 1,
    box_size::Unitful.Quantity = 1000UnitfulAstro.kpc,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    format::String = ".png",
)::Nothing

    # Get the simulation data.
    sim = getSnapshots(base_name, source_path)
    snap_files = sim["snap_files"]
    snap_numbers = sim["numbers"]

    # Create a directory to store the figures, if it doesn't exist.
    mkpath(output_path * "images/")

    # Generate and store the plots.
    short_snaps = @view snap_files[1:step:end]                       
    animation = @animate for (i, snapshot) in enumerate(short_snaps)

        positions = positionData(snapshot; sim_cosmo, box_size, length_unit)

        # Snashot number.
        number = snap_numbers[1 + step * (i - 1)]

        savefig(
            scatterGridPlot(positions),
            output_path * "images/" * base_name * "_" * number * format,
        )

    end

    # Make the GIF.
    gif(animation, output_path * anim_name * ".gif", fps = frame_rate)

    # Make the video.
    makeVideo(output_path * "images/", format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    densityMapPipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64; 
        <keyword arguments>
    )::Nothing

Save the results of the densityMapPlot function as one image per snapshot, 
and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base name of the snapshot files, 
  set in the GADGET variable SnapshotFileBase.
- `source_path::String`: Path to the directory containing the snapshot files, 
  set in the GADGET variable OutputDir.
- `anim_name::String`: Name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `output_path::String = "density_map/"`: Path to the output directory. The images will 
  be stored in `output_path`images/ and will be named `base_name`_XXX`format` where XXX 
  is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `step::Int64 = 1`: Step used to traverse the list of snapshots. The default is 1, 
  i.e. all snapshots will be plotted.
- `plane::String = "All"`: Indicates which plane will be plotted. 
  "XY" -> XY plane alone.
  "XZ" -> XZ plane alone.
  "YZ" -> YZ plane alone.
  "All" -> The three planes in a single 1x3 figure.
- `axes::Bool = false`: If true, the axes pasing through (0.0, 0.0) are drawn. If false, 
  no axes are drawn.
- `box_size::Unitful.Quantity = 1000UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. It has to have units, e.g. 1000UnitfulAstro.kpc, 
  which is the default. Its units don't have to be the same as `length_unit`.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.kpc, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by GR 
  can be used, namely ".pdf", ".ps", ".svg" and ".png", which is the default. 
"""
function densityMapPipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64;
    output_path::String = "density_map/",
    sim_cosmo::Int64 = 0,
    step::Int64 = 1,
    plane::String = "All",
    axes::Bool = false,
    box_size::Unitful.Quantity = 1000UnitfulAstro.kpc,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    format::String = ".png",
)::Nothing

    # Get the simulation data.
    sim = getSnapshots(base_name, source_path)
    snap_files = sim["snap_files"]
    snap_numbers = sim["numbers"]

    # Create a directory to store the figures, if it doesn't exist.
    mkpath(output_path * "images/")

    # Generate and store the plots.
    short_snaps = @view snap_files[1:step:end]                       
    animation = @animate for (i, snapshot) in enumerate(short_snaps)

        pos = positionData(snapshot; sim_cosmo, box_size, length_unit)
        mass = massData(snapshot, "gas"; sim_cosmo)
        density = densityData(snapshot; sim_cosmo)
        hsml = hsmlData(snapshot; sim_cosmo)

        # Snashot number.
        number = snap_numbers[1 + step * (i - 1)]

        savefig(
            densityMapPlot(pos, mass, density, hsml; plane, axes),
            output_path * "images/" * base_name * "_" * number * format,
        )
        
    end

    # Make the GIF.
    gif(animation, output_path  * anim_name * ".gif", fps = frame_rate)

    # Make the video.
    makeVideo(output_path * "images/", format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    starMapPipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64; 
        <keyword arguments>
    )::Nothing

Save the results of the starMapPlot function as one image per snapshot, 
and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base name of the snapshot files, 
  set in the GADGET variable SnapshotFileBase.
- `source_path::String`: Path to the directory containing the snapshot files, 
  set in the GADGET variable OutputDir.
- `anim_name::String`: Name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `output_path::String = "star_map/"`: Path to the output directory. The images will 
  be stored in `output_path`images/ and will be named `base_name`_XXX`format` where XXX 
  is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `step::Int64 = 1`: Step used to traverse the list of snapshots. The default is 1, 
  i.e. all snapshots will be plotted.
- `plane::String = "All"`: Indicates which plane will be plotted. 
  "XY" -> XY plane alone.
  "XZ" -> XZ plane alone.
  "YZ" -> YZ plane alone.
  "All" -> The three planes in a single 1x3 figure.
- `box_size::Unitful.Quantity = 1000UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. It has to have units, e.g. 1000UnitfulAstro.kpc, 
  which is the default. Its units don't have to be the same as `length_unit`.
- `box_factor::Float64 = 1.0`: Multiplicative factor for the plotting region. 
  It will scale `positions["box_size"]` if vacuum boundary conditions were used, and
  it will scale `positions["box_size"] / 2` if periodic boundary conditions were used.
- `axes::Bool = false`: If true, the axes pasing through (0.0, 0.0) are drawn. If false, 
  no axes are drawn.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.kpc, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by GR 
  can be used, namely ".pdf", ".ps", ".svg" and ".png", which is the default. 
"""
function starMapPipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64;
    output_path::String = "star_map/",
    sim_cosmo::Int64 = 0,
    step::Int64 = 1,
    plane::String = "All",
    box_size::Unitful.Quantity = 1000UnitfulAstro.kpc,
    box_factor::Float64 = 1.0, 
    axes::Bool = false,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,  
    format::String = ".png",
)::Nothing

    # Get the simulation data.
    sim = getSnapshots(base_name, source_path)
    snap_files = sim["snap_files"]
    snap_numbers = sim["numbers"]

    # Create a directory to store the figures, if it doesn't exist.
    mkpath(output_path * "images/")

    # Generate and store the plots.
    short_snaps = @view snap_files[1:step:end]                       
    animation = @animate for (i, snapshot) in enumerate(short_snaps)

        pos = positionData( snapshot; sim_cosmo, box_size, length_unit)

        # Snashot number.
        number = snap_numbers[1 + step * (i - 1)]

        savefig(
            starMapPlot(pos; plane, box_factor, axes),
            output_path * "images/" * base_name * "_" * number * format,
        )

    end

    # Make the GIF.
    gif(animation, output_path * anim_name * ".gif", fps = frame_rate)

    # Make the video.
    makeVideo( output_path * "images/", format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    gasStarEvolutionPipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64; 
        <keyword arguments>
    )::Nothing

Save the results of the gasStarEvolutionPlot for the last snapshot as one image and 
generate a GIF and a video animating the whole evolution for all snapshots. 
                                
# Arguments
- `base_name::String`: Base name of the snapshot files, 
  set in the GADGET variable SnapshotFileBase.
- `source_path::String`: Path to the directory containing the snapshot files, 
  set in the GADGET variable OutputDir.
- `anim_name::String`: Name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `output_path::String = "gas_star_evolution/"`: Path to the output directory. The image 
  will be stored in `output_path`images/ and will be named `base_name`_XXX`format` where XXX 
  is the number of the last snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `step::Int64 = 1`: Step used to traverse the list of snapshots. The default is 1, 
  i.e. all snapshots will be plotted.
- `box_size::Unitful.Quantity = 1000UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. It has to have units, e.g. 1000UnitfulAstro.kpc, 
  which is the default. Its units don't have to be the same as `length_unit`.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.kpc, which is the default.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Myr, which is the default.
- `sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr`: Unit of mass/time to 
  be used in the output, all available time and mass units in Unitful and UnitfulAstro 
  can be used, e.g. UnitfulAstro.Msun/UnitfulAstro.yr, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by GR 
  can be used, namely ".pdf", ".ps", ".svg" and ".png", which is the default. 
"""
function gasStarEvolutionPipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64;
    output_path::String = "gas_star_evolution/",
    sim_cosmo::Int64 = 0,
    step::Int64 = 1,
    box_size::Unitful.Quantity = 1000UnitfulAstro.kpc,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr,
    format::String = ".png",
)::Nothing

    # Get the simulation data.
    sim = getSnapshots(base_name, source_path)
    snap_files = sim["snap_files"]
    snap_numbers = sim["numbers"]

    # Create a directory to store the figures, if it doesn't exist.
    mkpath(output_path * "TEMP/")

    time_series = timeSeriesData(snap_files; sim_cosmo, time_unit, sfr_unit)

    # Generate and store the plots. 
    short_snaps = @view snap_files[1:step:end]                       
    animation = @animate for (i, snapshot) in enumerate(short_snaps)

        positions = positionData(snapshot; sim_cosmo, box_size, length_unit)

        # Snashot number.
        number = snap_numbers[1 + step * (i - 1)]

        savefig(
            gasStarEvolutionPlot(time_series, positions, i),
            output_path * "TEMP/" * base_name * "_" * number * format,
        )

    end

    # Make the GIF.
    gif(animation, output_path * anim_name * ".gif", fps = frame_rate)

    # Make the video.
    makeVideo(output_path * "TEMP/", format, output_path, anim_name, frame_rate)

    # Move the last figure out of the temporary directory.
    mv(
        output_path * "TEMP/" * base_name * "_" * snap_numbers[end] * format,
        output_path * anim_name * format,
        force = true,
    )

    # Delete the temporary directory and all its contents.
    rm(output_path * "TEMP/", recursive = true)

    return nothing
end

"""
    CMDFPipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64; 
        <keyword arguments>
    )::Nothing

Save the results of the CMDFPlot function as one image per snapshot 
(if there are stars present), and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable SnapshotFileBase.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable OutputDir.
- `anim_name::String`: Name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `output_path::String = "CMDF/"`: Path to the output directory. The images 
  will be stored in `output_path`images/ and will be named `base_name`_XXX`format` where XXX 
  is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `step::Int64 = 1`: Step used to traverse the list of snapshots. The default is 1, 
  i.e. all snapshots will be plotted.
- `x_norm::Bool = false`: If the x axis will be normalize to its maximum value. 
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Myr, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by 
  pgfplotsx can be used, namely ".pdf", ".tex", ".svg" and ".png", which is the default. 
"""
function CMDFPipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64;
    output_path::String = "CMDF/",
    sim_cosmo::Int64 = 0,
    step::Int64 = 1,
    x_norm::Bool = false,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    format::String = ".png",
)::Nothing

    # Create a directory to store the figures, if it doesn't exist.
    mkpath(output_path * "images/")

    # Get the simulation data.
    sim = getSnapshots(base_name, source_path)
    snap_files = sim["snap_files"]
    snap_numbers = sim["numbers"]

    time_data = timeSeriesData(snap_files; sim_cosmo)

    # Generate and store the plots.
    short_snaps = @view snap_files[1:step:end] 
    # animation = @animate          
    for (i, snapshot) in enumerate(short_snaps)

        header = read_header(snapshot)
        
        if header.nall[5] != 0
            mass_data = massData(snapshot, "stars"; sim_cosmo)
            z_data = zData(snapshot, "stars"; sim_cosmo)

            # Snashot number.
            number = snap_numbers[1 + step * (i - 1)]

            savefig(
                CMDFPlot(
                    mass_data, 
                    z_data,
                    time_data["clock_time"][1 + step * (i - 1)] * time_unit; 
                    bins = 50, 
                    x_norm
                ),
                output_path * "images/" * base_name * "_" * number * format,
            )
        end

    end

    # Make the GIF.
    # gif(animation, output_path * anim_name * ".gif", fps = frame_rate)

    # Make the video.
    # makeVideo(output_path * "images/", format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    CMDFPipeline(
        base_name::Vector{String},
        source_path::Vector{String},
        anim_name::String,
        frame_rate::Int64,
        labels::Array{String, 2}; 
        <keyword arguments>
    )::Nothing

Save the results of the CMDFPlot function for several simulations as one image per snapshot 
(if there are stars present), and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::Vector{String}`: Base names of the snapshot files, set in the GADGET 
  variable SnapshotFileBase.
- `source_path::Vector{String}`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable OutputDir.
- `anim_name::String`: Name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `labels::Array{String, 2}`: Labels for the different simulations.
- `output_path::String = "CMDF/"`: Path to the output directory. The images 
  will be stored in `output_path`images/ and will be named `base_name`_XXX`format` where XXX 
  is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `step::Int64 = 1`: Step used to traverse the list of snapshots. The default is 1, 
  i.e. all snapshots will be plotted.
- `x_norm::Bool = false`: If the x axis will be normalize to its maximum value. 
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Myr, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by 
  pgfplotsx can be used, namely ".pdf", ".tex", ".svg" and ".png", which is the default. 
"""
function CMDFPipeline(
    base_name::Vector{String},
    source_path::Vector{String},
    anim_name::String,
    frame_rate::Int64,
    labels::Array{String, 2};
    output_path::String = "mass_profile/",
    sim_cosmo::Int64 = 0,
    step::Int64 = 1,
    x_norm::Bool = false,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    format::String = ".png",
)::Nothing

    # Create a directory to store the figures, if it doesn't exist.
    mkpath(output_path * "images/")

    # Get the simulation data.
    snap_files = [getSnapshots(base_name[i], path)["snap_files"] 
                for (i, path) in enumerate(source_path)]

    # Length of the shortest simulation.
    min_len = minimum(length.(snap_files))

    # Time stamps, it should be the same for every dataset.
    time_data = timeSeriesData(snap_files[1]; sim_cosmo, time_unit)

    # Generate and store the plots.
    # animation = @animate 
    for i in 1:step:min_len
        
        headers = [read_header(snapshots[i]) for snapshots in snap_files]
        num_stars = getindex.(getfield.(headers, :nall), 5)

        if all(num_stars .!= 0)
    
            masses = [
                massData(snapshots[i], "stars"; sim_cosmo) 
                for snapshots in snap_files
            ]
            metallicities = [
                zData(snapshots[i], "stars"; sim_cosmo) 
                for snapshots in snap_files
            ]

            savefig(
                CMDFPlot(
                    masses, 
                    metallicities,
                    time_data["clock_time"][i] * time_unit,
                    labels;
                    bins = 50, 
                    x_norm
                ),
                output_path * "images/frame_" * string(i - 1) * format,
            )

        end

    end

    # Make the GIF.
    # gif(animation, output_path * anim_name * ".gif", fps = frame_rate)

    # Make the video.
    # makeVideo(output_path * "images/", format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    birthHistogramPipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64; 
        <keyword arguments>
    )::Nothing

Save the results of the birthHistogramPlot function as one image per snapshot 
(if there are stars present), and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable SnapshotFileBase.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable OutputDir.
- `anim_name::String`: Name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `output_path::String = "birth_histogram/"`: Path to the output directory. The images 
  will be stored in `output_path`images/ and will be named `base_name`_XXX`format` where XXX 
  is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `step::Int64 = 1`: Step used to traverse the list of snapshots. The default is 1, 
  i.e. all snapshots will be plotted.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.kpc, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by 
  pgfplotsx can be used, namely ".pdf", ".tex", ".svg" and ".png", which is the default. 
"""
function birthHistogramPipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64;
    output_path::String = "birth_histogram/",
    sim_cosmo::Int64 = 0,
    step::Int64 = 1,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    format::String = ".png",
)::Nothing

    # Get the simulation data.
    sim = getSnapshots(base_name, source_path)
    snap_files = sim["snap_files"]
    snap_numbers = sim["numbers"]

    time_data = timeSeriesData(snap_files; sim_cosmo)

    # Create a directory to store the figures, if it doesn't exist.
    mkpath(output_path * "images/")

    # Generate and store the plots.
    short_snaps = @view snap_files[1:step:end] 
    
    # animation = @animate          
    for (i, snapshot) in enumerate(short_snaps)

        header = read_header(snapshot)

        if header.nall[5] != 0
            nursery = birthPlace(
                i, 
                snap_files, 
                time_data["clock_time"],
                time_data["units"]["time"];
                sim_cosmo, 
                length_unit, 
                time_unit = time_data["units"]["time"],
            )

            # Snashot number.
            number = snap_numbers[1 + step * (i - 1)]

            savefig(
                birthHistogramPlot(nursery, bins = 50),
                output_path * "images/" * base_name * "_" * number * format,
            )
        end
    end

    # Make the GIF.
    # gif(animation, output_path * anim_name * ".gif", fps = frame_rate)

    # Make the video.
    # makeVideo(output_path * "images/", format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    evolutionSummaryPipeline(
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
  set in the GADGET variable SnapshotFileBase.
- `source_path::String`: Path to the directory containing the snapshot files, 
  set in the GADGET variable OutputDir.
- `fig_name::String`: Base name for the figures. The images will be named
  `fig_name`_vs_XXX`format` where XXX is 'time', 'redshift' or 'scale_factor'.
- `output_path::String = "evolution_summary/"`: Path to the output directory. The images 
  will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `mass_factor::Int64 = 0`: Numerical exponent to scale the mass, e.g. if mass_factor = 10 
  the corresponding axis will be scaled by 10^10.
- `number_factor::Int64 = 0`: Numerical exponent to scale the number of particles, 
  e.g. if number_factor = 4 the corresponding axis will be scaled by 10^4.
- `mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun`: Unit of mass to be used in the output, 
  all available mass units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Msun, which is the default.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Myr, which is the default.
- `sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr`: Unit of mass/time to 
  be used in the output, all available time and mass units in Unitful and UnitfulAstro 
  can be used, e.g. UnitfulAstro.Msun/UnitfulAstro.yr, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by 
  pgfplotsx can be used, namely ".pdf", ".tex", ".svg" and ".png", which is the default. 
"""
function evolutionSummaryPipeline(
    base_name::String,
    source_path::String,
    fig_name::String;
    output_path::String = "evolution_summary/",
    sim_cosmo::Int64 = 0,
    mass_factor::Int64 = 0,
    number_factor::Int64 = 0,
    mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr,
    format::String = ".png",
)::Nothing

    snap_files = getSnapshots(base_name, source_path)

    # Create a directory to store the figures, if it doesn't exist.
    mkpath(output_path)

    time_series = timeSeriesData(
                    snap_files["snap_files"]; 
                    sim_cosmo, 
                    mass_unit, 
                    time_unit, 
                    sfr_unit,
                )

    # Parameters vs. time. 
    savefig(
        timeSeriesPlot(time_series; mass_factor, number_factor),
        output_path * fig_name * "_vs_time" * format,
    )

    if sim_cosmo == 1

        # Parameters vs. scale factor. 
        savefig(
            scaleFactorSeriesPlot(time_series;mass_factor, number_factor),
            output_path * fig_name * "_vs_scale_factor" * format,
        )
        # Parameters vs. redshift. 
        savefig(
            redshiftSeriesPlot(time_series; mass_factor, number_factor),
            output_path * fig_name * "_vs_redshift" * format,
        )

    end

    return nothing
end

"""
    compareSimulationsPipeline(
        base_name::Vector{String},
        source_path::Vector{String},
        labels::Array{String, 2},
        fig_name::String,
        x_quantity::String,
        y_quantity::String; 
        <keyword arguments>
    )::Nothing

Make a figure comparing `y_quantity` vs. `x_quantity` for several simulations.

`x_quantity` and `y_quantity` can be any magnitude used in the timeSeriesData 
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
- "gas_bar_frac" (Baryonic gas fraction)                  
- "star_bar_frac" (Baryonic star fraction)

# Arguments
- `base_name::Vector{String}`: Base names of the snapshot files, set in the GADGET 
  variable SnapshotFileBase.
- `source_path::Vector{String}`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable OutputDir.
- `labels::Array{String, 2}`: Labels for the different simulations, e.g. [label1 label2 ...].
- `fig_name::String`: Base name for the figure. The file will be named
  `fig_name`_`y_quantity`_vs_`x_quantity` `format`.
- `x_quantity::String`: Physical magnitude for the x axis. 
- `y_quantity::String`: Physical magnitude for the y axis.
- `output_path::String = "compare_simulations/"`: Path to the output directory. The images 
  will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `title::String = ""`: Title for the figure. If an empty string is given no title is 
  printed, which is the default.
- `x_factor::Int64 = 0`: Numerical exponent to scale the `x_quantity`, e.g. if x_factor = 10 
  the corresponding axis will be scaled by 10^10. The default is 0, i.e. no scaling.
- `y_factor::Int64 = 0`: Numerical exponent to scale the `y_quantity`, e.g. if y_factor = 10 
  the corresponding axis will be scaled by 10^10. The default is 0, i.e. no scaling.
- `scale::NTuple{2, Symbol} = (:identity, :identity)`: Scaling to be used for the x 
  and y axes. The two options are:
  :identity => no scaling.
  :log10 => logarithmic scaling.
- `smooth_data::Bool = false`: If true a smoothing window with no weighs is applied to 
  the y data. If false (the default) no transformation occurs.
- `bins::Int64 = 0`: Number of subdivisions for the smoothing of the data, only relevant if
  `smooth_data = true`. 
- `legend_pos::Symbol = :bottomright`: Position of the legend, e.g. :topleft.
- `mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun`: Unit of mass to be used in the output, 
  all available mass units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Msun, which is the default.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Myr, which is the default.
- `sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr`: Unit of mass/time to 
  be used in the output, all available time and mass units in Unitful and UnitfulAstro 
  can be used, e.g. UnitfulAstro.Msun/UnitfulAstro.yr, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by 
  pgfplotsx can be used, namely ".pdf", ".tex", ".svg" and ".png", which is the default. 
"""
function compareSimulationsPipeline(
    base_name::Vector{String},
    source_path::Vector{String},
    labels::Array{String, 2},
    fig_name::String,
    x_quantity::String,
    y_quantity::String;
    output_path::String = "compare_simulations/",
    sim_cosmo::Int64 = 0,
    title::String = "",
    x_factor::Int64 = 0,
    y_factor::Int64 = 0,
    scale::NTuple{2, Symbol} = (:identity, :identity),
    smooth_data::Bool = false, 
    bins::Int64 = 50,
    legend_pos::Symbol = :bottomright,
    mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr,
    format::String = ".png",
)::Nothing

    time_series = [timeSeriesData(
                        getSnapshots(name, path)["snap_files"]; 
                        sim_cosmo, 
                        mass_unit, 
                        time_unit, 
                        sfr_unit,
                    ) for (name, path) in zip(base_name, source_path)]

    # Create a directory to store the figure, if it doesn't exist.
    mkpath(output_path)

    # `y_quantity` vs. `x_quantity` plot.  
    savefig(
        compareSimulationsPlot(
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
        ), 
        output_path * fig_name * "_" * y_quantity * "_vs_" * x_quantity * format,
    )

    return nothing
end

"""
    densityHistogramPipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64; 
        <keyword arguments>
    )::Nothing

Save the results of the densityHistogramPlot function as one image per snapshot, 
and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable SnapshotFileBase.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable OutputDir.
- `anim_name::String`: Name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `output_path::String = "density_histogram/"`: Path to the output directory. The images 
  will be stored in `output_path`images/ and will be named `base_name`_XXX`format` where XXX 
  is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `step::Int64 = 1`: Step used to traverse the list of snapshots. The default is 1, 
  i.e. all snapshots will be plotted.
- `factor::Int64 = 0`: Numerical exponent to scale the density, e.g. if factor = 10 
  the y axis will be scaled by 10^10. The default is 0, i.e. no scaling.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Myr, which is the default.
- `density_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.kpc^3`: Unit of density 
  to be used in the output, all available density units in Unitful and UnitfulAstro can 
  be used, e.g. UnitfulAstro.Msun / UnitfulAstro.kpc^3, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by 
  pgfplotsx can be used, namely ".pdf", ".tex", ".svg" and ".png", which is the default. 
"""
function densityHistogramPipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64;
    output_path::String = "density_histogram/",
    sim_cosmo::Int64 = 0,
    step::Int64 = 1,
    factor::Int64 = 0,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    density_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.kpc^3,
    format::String = ".png",
)::Nothing

    # Get the simulation data.
    sim = getSnapshots(base_name, source_path)
    snap_files = sim["snap_files"]
    snap_numbers = sim["numbers"]

    time_data = timeSeriesData(snap_files; sim_cosmo, time_unit)

    # Create a directory to store the figures, if it doesn't exist.
    mkpath(output_path * "images/")

    # Generate and store the plots.
    short_snaps = @view snap_files[1:step:end]
    # animation = @animate 
    for (i, snapshot) in enumerate(short_snaps)

        density = densityData(snapshot; sim_cosmo)

        # Snashot number.
        number = snap_numbers[1 + step * (i - 1)]

        savefig(
            densityHistogramPlot(
                density,
                time_data["clock_time"][1 + step * (i - 1)] * time_unit;
                factor,
            ), 
            output_path * "images/" * base_name * "_" * number * format,
        )

    end

    # Make the GIF.
    # gif(animation, output_path  * anim_name * ".gif", fps = frame_rate)

    # Make the video.
    # makeVideo(output_path * "images/", format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    densityProfilePipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64,
        type::String; 
        <keyword arguments>)::Nothing

Save the results of the densityProfilePlot function as one image per snapshot,
and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable SnapshotFileBase.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable OutputDir.
- `anim_name::String`: Name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `type::String`: Particle type.
  "gas" -> Gas particle. 
  "dark_matter" -> Dark matter particle.
  "stars" -> Star particle.
- `output_path::String = "density_profile/"`: Path to the output directory. The images 
  will be stored in `output_path`images/ and will be named `base_name`_XXX`format` where XXX 
  is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `scale::Symbol = :identity`: Scaling to be used for the y axis.
  The two options are:
  :identity => no scaling.
  :log10 => logarithmic scaling.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. The default is 1, 
  i.e. all snapshots will be plotted.
- `bins::Int64 = 100`: Number of subdivisions of the region to be used for the profile. 
  The default is 100.
- `factor::Int64 = 0`: Numerical exponent to scale the density, e.g. if factor = 10 
  the y axis will be scaled by 10^10. The default is 0, i.e. no scaling.
- `box_factor::Float64 = 1.0`: Multiplicative factor for the plotting region. 
  It will scale `positions["box_size"]` if vacuum boundary conditions were used, and
  it will scale `positions["box_size"] / 2` if periodic boundary conditions were used.
- `box_size::Unitful.Quantity = 1000UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. It has to have units, e.g. 1000UnitfulAstro.kpc, 
  which is the default. Its units don't have to be the same as `length_unit`.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.kpc, which is the default.
- `density_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.kpc^3`: Unit of density 
  to be used in the output, all available density units in Unitful and UnitfulAstro can 
  be used, e.g. UnitfulAstro.Msun / UnitfulAstro.kpc^3, which is the default.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Myr, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by 
  pgfplotsx can be used, namely ".pdf", ".tex", ".svg" and ".png", which is the default. 
"""
function densityProfilePipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64,
    type::String;
    output_path::String = "density_profile/",
    sim_cosmo::Int64 = 0,
    scale::Symbol = :identity,
    step::Int64 = 1,
    bins::Int64 = 100,
    factor::Int64 = 0,
    box_factor::Float64 = 1.0,
    box_size::Unitful.Quantity = 1000UnitfulAstro.kpc,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    density_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.kpc^3,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    format::String = ".png",
)::Nothing

    # Get the simulation data.
    sim = getSnapshots(base_name, source_path)
    snap_files = sim["snap_files"]
    snap_numbers = sim["numbers"]

    time_data = timeSeriesData(snap_files; sim_cosmo, time_unit)

    # Create a directory to store the figures, if it doesn't exist.
    mkpath(output_path * "images/")

    # Generate and store the plots.
    short_snaps = @view snap_files[1:step:end] 
    # animation = @animate          
    for (i, snapshot) in enumerate(short_snaps)

        positions = positionData(snapshot; sim_cosmo, box_size, length_unit)
        mass = massData(snapshot, type; sim_cosmo)

        # Snashot number.
        number = snap_numbers[1 + step * (i - 1)]

        savefig(
            densityProfilePlot(
                positions,
                mass,
                time_data["clock_time"][1 + step * (i - 1)] * time_unit;
                scale,
                bins,
                factor,
                box_factor,
            ),
            output_path * "images/" * base_name * "_" * number * format,
        )

    end

    # Make the GIF.
    # gif(animation, output_path * anim_name * ".gif", fps = frame_rate)

    # Make the video.
    # makeVideo(output_path * "images/", format, output_path * type, anim_name, frame_rate)

    return nothing
end

"""
    densityProfilePipeline(
        base_name::Vector{String},
        source_path::Vector{String},
        anim_name::String,
        frame_rate::Int64,
        type::String,
        labels::Array{String, 2}; 
        <keyword arguments>
    )::Nothing

Save the results of the densityProfilePlot function for several simulations as one image 
per snapshot, and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::Vector{String}`: Base names of the snapshot files, set in the GADGET 
  variable SnapshotFileBase.
- `source_path::Vector{String}`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable OutputDir.
- `anim_name::String`: Name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `type::String`: Particle type.
  "gas" -> Gas particle. 
  "dark_matter" -> Dark matter particle.
  "stars" -> Star particle.
- `labels::Array{String, 2}`: Labels for the different simulations.
- `output_path::String = "density_profile/"`: Path to the output directory. The images 
  will be stored in `output_path`images/ and will be named frame_XXX`format` where XXX is 
  the ordinal of the frame. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `scale::Symbol = :identity`: Scaling to be used for the y axis.
  The two options are:
  :identity => no scaling.
  :log10 => logarithmic scaling.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. The default is 1, 
  i.e. all snapshots will be plotted.
- `bins::Int64 = 100`: Number of subdivisions of the region to be used for the profile. 
  The default is 100.
- `factor::Int64 = 0`: Numerical exponent to scale the density, e.g. if factor = 10 
  the y axis will be scaled by 10^10. The default is 0, i.e. no scaling.
- `box_factor::Float64 = 1.0`: Multiplicative factor for the plotting region. 
  It will scale `positions["box_size"]` if vacuum boundary conditions were used, and
  it will scale `positions["box_size"] / 2` if periodic boundary conditions were used.
- `box_size::Unitful.Quantity = 1000UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. It has to have units, e.g. 1000UnitfulAstro.kpc, 
  which is the default. Its units don't have to be the same as `length_unit`.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.kpc, which is the default.
- `density_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.kpc^3`: Unit of density 
  to be used in the output, all available density units in Unitful and UnitfulAstro can 
  be used, e.g. UnitfulAstro.Msun / UnitfulAstro.kpc^3, which is the default.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Myr, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by 
  pgfplotsx can be used, namely ".pdf", ".tex", ".svg" and ".png", which is the default. 
"""
function densityProfilePipeline(
    base_name::Vector{String},
    source_path::Vector{String},
    anim_name::String,
    frame_rate::Int64,
    type::String,
    labels::Array{String, 2};
    output_path::String = "density_profile/",
    sim_cosmo::Int64 = 0,
    scale::Symbol = :identity,
    step::Int64 = 1,
    bins::Int64 = 100,
    factor::Int64 = 0,
    box_factor::Float64 = 1.0,
    box_size::Unitful.Quantity = 1000UnitfulAstro.kpc,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    density_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.kpc^3,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    format::String = ".png",
)::Nothing

    # Create a directory to store the figures, if it doesn't exist.
    mkpath(output_path * "images/")

    # Get the simulation data.
    snap_files = [getSnapshots(base_name[i], path)["snap_files"] 
                for (i, path) in enumerate(source_path)]
    
    # Length of the shortest simulation.
    min_len = minimum(length.(snap_files))

    # Time stamps, it should be the same for every dataset.
    time_data = timeSeriesData(snap_files[1]; sim_cosmo, time_unit)

    # Generate and store the plots.
    # animation = @animate 
    for i in 1:step:min_len

        positions = [positionData(snapshots[i]; sim_cosmo, box_size, length_unit) 
                    for snapshots in snap_files]
        masses = [massData(snapshots[i], type; sim_cosmo) for snapshots in snap_files]

        savefig(
            densityProfilePlot(
                positions,
                masses,
                time_data["clock_time"][i] * time_unit,
                labels;
                scale,
                bins,
                factor,
                box_factor,
            ),
            output_path * "images/frame_"  * string(i - 1) * format,
        )
    end

    # Make the GIF.
    # gif(animation, output_path * anim_name * ".gif", fps = frame_rate)

    # Make the video.
    # makeVideo(output_path * "images/", format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    metallicityProfilePipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64,
        type::String; 
        <keyword arguments>
    )::Nothing

Save the results of the metallicityProfilePlot function as one image per snapshot,
and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable SnapshotFileBase.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable OutputDir.
- `anim_name::String`: Name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `type::String`: Particle type.
  "gas" -> Gas particle. 
  "dark_matter" -> Dark matter particle.
  "stars" -> Star particle.
- `output_path::String = "metallicity_profile/"`: Path to the output directory. The images 
  will be stored in `output_path`images/ and will be named `base_name`_XXX`format` where XXX 
  is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `step::Int64 = 1`: Step used to traverse the list of snapshots. The default is 1, 
  i.e. all snapshots will be plotted.
- `bins::Int64 = 100`: Number of subdivisions of the region to be used for the profile. 
  The default is 100.
- `box_factor::Float64 = 1.0`: Multiplicative factor for the plotting region. 
  It will scale `positions["box_size"]` if vacuum boundary conditions were used, and
  it will scale `positions["box_size"] / 2` if periodic boundary conditions were used.
- `box_size::Unitful.Quantity = 1000UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. It has to have units, e.g. 1000UnitfulAstro.kpc, 
  which is the default. Its units don't have to be the same as `length_unit`.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.kpc, which is the default.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Myr, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by 
  pgfplotsx can be used, namely ".pdf", ".tex", ".svg" and ".png", which is the default. 
"""
function metallicityProfilePipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64,
    type::String;
    output_path::String = "metallicity_profile/",
    sim_cosmo::Int64 = 0,
    scale::Symbol = :identity,
    step::Int64 = 1,
    bins::Int64 = 100,
    box_factor::Float64 = 1.0,
    box_size::Unitful.Quantity = 1000UnitfulAstro.kpc,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    format::String = ".png",
)::Nothing

    # Create a directory to store the figures, if it doesn't exist.
    mkpath(output_path * "images/")

    # Get the simulation data.
    sim = getSnapshots(base_name, source_path)
    snap_files = sim["snap_files"]
    snap_numbers = sim["numbers"]

    time_data = timeSeriesData(snap_files; sim_cosmo, time_unit)

    # Generate and store the plots.
    short_snaps = @view snap_files[1:step:end] 
    # animation = @animate          
    for (i, snapshot) in enumerate(short_snaps)

        positions = positionData(snapshot; sim_cosmo, box_size, length_unit)
        mass = massData(snapshot, type; sim_cosmo)
        metallicities = zData(snapshot, type; sim_cosmo)

        # Snashot number.
        number = snap_numbers[1 + step * (i - 1)]

        savefig(
            metallicityProfilePlot(
                positions,
                mass,
                metallicities,
                time_data["clock_time"][1 + step * (i - 1)] * time_unit;
                scale,
                bins,
                box_factor,
            ),
            output_path * "images/" * base_name * "_" * number * format,
        )
    end

    # Make the GIF.
    # gif(animation, output_path * anim_name * ".gif", fps = frame_rate)

    # Make the video.
    # makeVideo(output_path * "images/", format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    metallicityProfilePipeline(
        base_name::Vector{String},
        source_path::Vector{String},
        anim_name::String,
        frame_rate::Int64,
        type::String,
        labels::Array{String, 2}; 
        <keyword arguments>
    )::Nothing

Save the results of the metallicityProfilePlot function for several simulations as one 
image per snapshot, and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::Vector{String}`: Base names of the snapshot files, set in the GADGET 
  variable SnapshotFileBase.
- `source_path::Vector{String}`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable OutputDir.
- `anim_name::String`: Name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `type::String`: Particle type.
  "gas" -> Gas particle. 
  "dark_matter" -> Dark matter particle.
  "stars" -> Star particle.
- `labels::Array{String,2}`: Labels for the different simulations.
- `output_path::String = "metallicity_profile/"`: Path to the output directory. The images 
  will be stored in `output_path`images/ and will be named frame_XXX`format` where XXX is 
  the ordinal of the frame. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `scale::Symbol = :identity`: Scaling to be used for the y axis.
  The two options are:
  :identity => no scaling.
  :log10 => logarithmic scaling.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. The default is 1, 
  i.e. all snapshots will be plotted.
- `bins::Int64=100`: Number of subdivisions of the region to be used for the profile. 
  The default is 100.
- `box_factor::Float64 = 1.0`: Multiplicative factor for the plotting region. 
  It will scale `positions["box_size"]` if vacuum boundary conditions were used, and
  it will scale `positions["box_size"] / 2` if periodic boundary conditions were used.
- `box_size::Unitful.Quantity = 1000UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. It has to have units, e.g. 1000UnitfulAstro.kpc, 
  which is the default. Its units don't have to be the same as `length_unit`.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.kpc, which is the default.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Myr, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by 
  pgfplotsx can be used, namely ".pdf", ".tex", ".svg" and ".png", which is the default. 
"""
function metallicityProfilePipeline(
    base_name::Vector{String},
    source_path::Vector{String},
    anim_name::String,
    frame_rate::Int64,
    type::String,
    labels::Array{String, 2};
    output_path::String = "metallicity_profile/",
    sim_cosmo::Int64 = 0,
    scale::Symbol = :identity,
    step::Int64 = 1,
    bins::Int64 = 100,
    box_factor::Float64 = 1.0,
    box_size::Unitful.Quantity = 1000UnitfulAstro.kpc,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,  
    format::String = ".png",
)::Nothing

    # Create a directory to store the figures, if it doesn't exist.
    mkpath(output_path * "images/")

    # Get the simulation data.
    snap_files = [getSnapshots(base_name[i], path)["snap_files"] 
                for (i, path) in enumerate(source_path)]

    # Length of the shortest simulation.
    min_len = minimum(length.(snap_files))

    # Time stamps, it should be the same for every dataset.
    time_data = timeSeriesData(snap_files[1]; sim_cosmo, time_unit)

    # Generate and store the plots.
    # animation = @animate 
    for i in 1:step:min_len

        positions = [positionData(snapshots[i]; sim_cosmo, box_size, length_unit) 
                    for snapshots in snap_files]
        masses = [massData(snapshots[i], type; sim_cosmo) for snapshots in snap_files]
        metallicities = [zData(snapshots[i], type; sim_cosmo) for snapshots in snap_files]

        savefig(
            metallicityProfilePlot(
                positions,
                masses,
                metallicities,
                time_data["clock_time"][i] * time_unit,
                labels;
                scale,
                bins,
                box_factor,
            ),
            output_path * "images/frame_" * string(i - 1) * format,
        )

    end

    # Make the GIF.
    # gif(animation, output_path * anim_name * ".gif", fps = frame_rate)

    # Make the video.
    # makeVideo(output_path * "images/", format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    massProfilePipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64,
        type::String; 
        <keyword arguments>
    )::Nothing

Save the results of the massProfilePlot function as one image per snapshot,
and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable SnapshotFileBase.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable OutputDir.
- `anim_name::String`: Name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `type::String`: Particle type.
  "gas" -> Gas particle. 
  "dark_matter" -> Dark matter particle.
  "stars" -> Star particle.
- `scale::Symbol = :identity`: Scaling to be used for the y axis.
  The two options are:
  :identity => no scaling.
  :log10 => logarithmic scaling.
- `output_path::String = "mass_profile/"`: Path to the output directory. The images 
  will be stored in `output_path`images/ and will be named `base_name`_XXX`format` where XXX 
  is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `step::Int64 = 1`: Step used to traverse the list of snapshots. The default is 1, 
  i.e. all snapshots will be plotted.
- `bins::Int64 = 100`: Number of subdivisions of the region to be used for the profile. 
  The default is 100.
- `factor::Int64 = 0`: Numerical exponent to scale the mass, e.g. if factor = 10 
  the y axis will be scaled by 10^10. The default is 0, i.e. no scaling.
- `box_factor::Float64 = 1.0`: Multiplicative factor for the plotting region. 
  It will scale `positions["box_size"]` if vacuum boundary conditions were used, and
  it will scale `positions["box_size"] / 2` if periodic boundary conditions were used.
- `box_size::Unitful.Quantity = 1000UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. It has to have units, e.g. 1000UnitfulAstro.kpc, 
  which is the default. Its units don't have to be the same as `length_unit`.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.kpc, which is the default.
- `mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun`: Unit of mass to be used in the output, 
  all available mass units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Msun, which is the default.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Myr, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by 
  pgfplotsx can be used, namely ".pdf", ".tex", ".svg" and ".png", which is the default. 
"""
function massProfilePipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64,
    type::String;
    output_path::String = "mass_profile/",
    sim_cosmo::Int64 = 0,
    scale::Symbol = :identity,
    step::Int64 = 1,
    bins::Int64 = 100,
    factor::Int64 = 0,
    box_factor::Float64 = 1.0,
    box_size::Unitful.Quantity = 1000UnitfulAstro.kpc,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    format::String = ".png",
)::Nothing

    # Get the simulation data.
    sim = getSnapshots(base_name, source_path)
    snap_files = sim["snap_files"]
    snap_numbers = sim["numbers"]

    time_data = timeSeriesData(snap_files; sim_cosmo, time_unit)

    # Create a directory to store the figures, if it doesn't exist.
    mkpath(output_path * "images/")

    # Generate and store the plots.
    short_snaps = @view snap_files[1:step:end] 
    # animation = @animate          
    for (i, snapshot) in enumerate(short_snaps)

        positions = positionData(snapshot; sim_cosmo, box_size, length_unit)
        mass = massData(snapshot, type; sim_cosmo)

        # Snashot number.
        number = snap_numbers[1 + step * (i - 1)]

        savefig(
            massProfilePlot(
                positions,
                mass,
                time_data["clock_time"][1 + step * (i - 1)] * time_unit;
                scale,
                bins,
                factor,
                box_factor,
            ),
            output_path * "images/" * base_name * "_" * number * format,
        )

    end

    # Make the GIF.
    # gif(animation, output_path * anim_name * ".gif", fps = frame_rate)

    # Make the video.
    # makeVideo(output_path * "images/", format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    massProfilePipeline(
        base_name::Vector{String},
        source_path::Vector{String},
        anim_name::String,
        frame_rate::Int64,
        type::String,
        labels::Array{String, 2}; 
        <keyword arguments>
    )::Nothing

Save the results of the massProfilePlot function for several simulations as one image 
per snapshot, and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::Vector{String}`: Base names of the snapshot files, set in the GADGET 
  variable SnapshotFileBase.
- `source_path::Vector{String}`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable OutputDir.
- `anim_name::String`: Name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `type::String`: Particle type.
  "gas" -> Gas particle. 
  "dark_matter" -> Dark matter particle.
  "stars" -> Star particle.
- `labels::Array{String, 2}`: Labels for the different simulations.
- `output_path::String = "mass_profile/"`: Path to the output directory. The images 
  will be stored in `output_path`images/ and will be named frame_XXX`format` where XXX is 
  the ordinal of the frame. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `scale::Symbol = :identity`: Scaling to be used for the y axis.
  The two options are:
  :identity => no scaling.
  :log10 => logarithmic scaling.
- `step::Int64 = 1`: Step used to traverse the list of snapshots. The default is 1, 
  i.e. all snapshots will be plotted.
- `bins::Int64 = 100`: Number of subdivisions of the region to be used for the profile. 
  The default is 100.
- `factor::Int64 = 0`: Numerical exponent to scale the mass, e.g. if factor = 10 
  the y axis will be scaled by 10^10. The default is 0, i.e. no scaling.
- `box_factor::Float64 = 1.0`: Multiplicative factor for the plotting region. 
  It will scale `positions["box_size"]` if vacuum boundary conditions were used, and
  it will scale `positions["box_size"] / 2` if periodic boundary conditions were used.
- `box_size::Unitful.Quantity = 1000UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. It has to have units, e.g. 1000UnitfulAstro.kpc, 
  which is the default. Its units don't have to be the same as `length_unit`.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.kpc, which is the default.
- `mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun`: Unit of mass to be used in the output, 
  all available mass units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Msun, which is the default.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Myr, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by 
  pgfplotsx can be used, namely ".pdf", ".tex", ".svg" and ".png", which is the default. 
"""
function massProfilePipeline(
    base_name::Vector{String},
    source_path::Vector{String},
    anim_name::String,
    frame_rate::Int64,
    type::String,
    labels::Array{String, 2};
    output_path::String = "mass_profile/",
    sim_cosmo::Int64 = 0,
    scale::Symbol = :identity,
    step::Int64 = 1,
    bins::Int64 = 100,
    factor::Int64 = 0,
    box_factor::Float64 = 1.0,
    box_size::Unitful.Quantity = 1000UnitfulAstro.kpc,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    format::String = ".png",
)::Nothing

    # Create a directory to store the figures, if it doesn't exist.
    mkpath(output_path * "images/")

    # Get the simulation data.
    snap_files = [getSnapshots(base_name[i], path)["snap_files"] 
                for (i, path) in enumerate(source_path)]

    # Length of the shortest simulation.
    min_len = minimum(length.(snap_files))

    # Time stamps, it should be the same for every dataset.
    time_data = timeSeriesData(snap_files[1]; sim_cosmo, time_unit)

    # Generate and store the plots.
    # animation = @animate 
    for i in 1:step:min_len

        positions = [positionData(snapshots[i]; sim_cosmo, box_size, length_unit) 
                    for snapshots in snap_files]
        masses = [massData(snapshots[i], type; sim_cosmo) for snapshots in snap_files]

        savefig(
            massProfilePlot(
                positions,
                masses,
                time_data["clock_time"][i] * time_unit,
                labels;
                scale,
                bins,
                factor,
                box_factor,
            ),
            output_path * "images/frame_" * string(i - 1) * format,
        )

    end

    # Make the GIF.
    # gif(animation, output_path * anim_name * ".gif", fps = frame_rate)

    # Make the video.
    # makeVideo(output_path * "images/", format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    sfrTxtPipeline(
        base_names::Vector{String},
        source_paths::Vector{String},
        x_axis::Int64,
        y_axis::Vector{Int64}; 
        <keyword arguments>
    )::Plots.Plot

Save the results of the sfrTxtPlot function as one image per simulation.

# Warning 
THIS FUNCTION NEEDS A FILE (sfr.txt) WHICH IS NOT PRODUCED BY ANY PUBLIC VERSION OF GADGET.

# Arguments
- `base_name::Vector{String}`: Base names of the snapshot files, set in the GADGET 
  variable SnapshotFileBase.
- `source_path::Vector{String}`: Paths to the directories containing the sfr.txt files, 
  set in the GADGET variable OutputDir.
- `x_axis::Int64`: Column number for the x axis.
- `y_axis::Vector{Int64}`: Vector of columns numbers for the y axis.
- `output_path::String = "sfr_txt/"`: Path to the output directory.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `title::Vector{String} = String[]`: Titles for the figures. If an empty string is given
  no title is printed, which is the default.
- `names::Vector{String} = String[]`: Names for the files. If an empty string is given the images
  will be asign a number, given by the order of `source_path`.
- `bins::Int64 = 0`: Number of subdivisions for the smoothing of the data. 
  The default is 0, i.e. no smoothing. It will apply equaly to every figure produced.
- `scale::NTuple{2, Symbol} = (:identity, :identity)`: Scaling to be used for the x and y 
  axes. It will apply equaly to every figure produced.
  The two options are:
  :identity => no scaling.
  :log10 => logarithmic scaling.
- `min_filter::NTuple{2, Float64} = (-Inf, -Inf)`: Value filter for the x and y axes. 
  It will apply equaly to every figure produced. If a value of the x data is lower 
  than min_filter[1], then it is deleted. Equivalently with the y axis and min_filter[2]. 
  The default is -Inf for both, i.e. no filtering.
- `mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun`: Unit of mass to be used in the output, 
  all available mass units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Msun, which is the default.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Myr, which is the default.
- `sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr`: Unit of mass/time to 
  be used in the output, all available time and mass units in Unitful and UnitfulAstro 
  can be used, e.g. UnitfulAstro.Msun/UnitfulAstro.yr, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by 
  pgfplotsx can be used, namely ".pdf", ".tex", ".svg" and ".png", which is the default. 

# Returns
- The plot generated by the PGFPlotsX backend of Plots.jl.
"""
function sfrTxtPipeline(
    base_name::Vector{String},
    source_path::Vector{String},
    x_axis::Int64,
    y_axis::Vector{Int64};
    output_path::String = "sfr_txt/",
    sim_cosmo::Int64 = 0,
    title::Vector{String} = String[],
    names::Vector{String} = String[],
    bins::Int64 = 0,
    scale::NTuple{2, Symbol} = (:identity, :identity),
    min_filter::NTuple{2, Float64} = (-Inf, -Inf),
    mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr,
    format::String = ".png",
)::Nothing

    # Create a directory to store the figures, if it doesn't exist.
    mkpath(output_path)

    # By default every figure has no title.
    if isempty(title)
        title = ["" for _ in source_path]
    end

    # By default every figure has a number as its name.
    if isempty(names)
        names = [string(i - 1) for i in eachindex(source_path)]
    end

    @inbounds for i in eachindex(source_path, base_name)

        sfr_data = sfrTxtData(
            source_path[i], 
            base_name[i]; 
            sim_cosmo, 
            mass_unit, 
            time_unit, 
            sfr_unit
        )

        savefig(
                sfrTxtPlot(
                    sfr_data, 
                    x_axis, 
                    y_axis;
                    title = title[i], 
                    bins, 
                    scale, 
                    min_filter,
                ),
                output_path * names[i] * format,
            )

    end

    return nothing
end

"""
    temperatureHistogramPipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64; 
        <keyword arguments>
    )::Nothing

Save the results of the temperatureHistogramPlot function as one image per snapshot, 
and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable SnapshotFileBase.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable OutputDir.
- `anim_name::String`: Name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `output_path::String = "temperature_histogram/"`: Path to the output directory. The images 
  will be stored in `output_path`images/ and will be named `base_name`_XXX`format` where XXX 
  is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `step::Int64 = 1`: Step used to traverse the list of snapshots. The default is 1, 
  i.e. all snapshots will be plotted.
- `temp_unit::Unitful.FreeUnits = Unitful.K`: Unit of temperature to be used in the 
  output, all available temperature units in Unitful and UnitfulAstro can be used, 
  e.g. Unitful.K, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by 
  pgfplotsx can be used, namely ".pdf", ".tex", ".svg" and ".png", which is the default. 
"""
function temperatureHistogramPipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64;
    output_path::String = "temperature_histogram/",
    sim_cosmo::Int64 = 0,
    step::Int64 = 1,
    temp_unit::Unitful.FreeUnits = Unitful.K,
    format::String = ".png",
)::Nothing

    # Get the simulation data.
    sim = getSnapshots(base_name, source_path)
    snap_files = sim["snap_files"]
    snap_numbers = sim["numbers"]

    time_data = timeSeriesData(snap_files; sim_cosmo)
    time_unit = time_data["units"]["time"]
    clock = time_data["clock_time"]

    # Create a directory to store the figures, if it doesn't exist.
    mkpath(output_path * "images/")

    # Generate and store the plots.
    short_snaps = @view snap_files[1:step:end] 
    
    # animation = @animate          
    for (i, snapshot) in enumerate(short_snaps)

        temp_data = tempData(snapshot; sim_cosmo, temp_unit)

        # Snashot number.
        number = snap_numbers[1 + step * (i - 1)]

        savefig(
            temperatureHistogramPlot(
                temp_data, 
                clock[1 + step * (i - 1)] * time_unit, 
                bins = 30
            ),
            output_path * "images/" * base_name * "_" * number * format,
        )
        
    end

    # Make the GIF.
    # gif(animation, output_path * anim_name * ".gif", fps = frame_rate)

    # Make the video.
    # makeVideo(output_path * "images/", format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    rhoTempPipeline(
        base_name::String,
        source_path::String,
        anim_name::String,
        frame_rate::Int64; 
        <keyword arguments>
    )::Nothing

Save the results of the rhoTempPlot function as one image per snapshot, 
and then generate a GIF and a video animating the images. 

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable SnapshotFileBase.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable OutputDir.
- `anim_name::String`: Name of the generated video and GIF, without the extension.
- `frame_rate::Int64`: Frame rate of the output video and GIF.
- `output_path::String = "rho_vs_temp/"`: Path to the output directory. The images 
  will be stored in `output_path`images/ and will be named `base_name`_XXX`format` where XXX 
  is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `step::Int64 = 1`: Step used to traverse the list of snapshots. The default is 1, 
  i.e. all snapshots will be plotted.
- `temp_unit::Unitful.FreeUnits = Unitful.K`: Unit of temperature to be used in the 
  output, all available temperature units in Unitful and UnitfulAstro can be used, 
  e.g. Unitful.K, which is the default.
- `density_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.kpc^3`: Unit of density 
  to be used in the output, all available density units in Unitful and UnitfulAstro can 
  be used, e.g. UnitfulAstro.Msun / UnitfulAstro.kpc^3, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by 
  pgfplotsx can be used, namely ".pdf", ".tex", ".svg" and ".png", which is the default. 
"""
function rhoTempPipeline(
    base_name::String,
    source_path::String,
    anim_name::String,
    frame_rate::Int64;
    output_path::String = "rho_vs_temp/",
    sim_cosmo::Int64 = 0,
    step::Int64 = 1,
    temp_unit::Unitful.FreeUnits = Unitful.K,
    density_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.kpc^3,
    format::String = ".png",
)::Nothing

    # Get the simulation data.
    sim = getSnapshots(base_name, source_path)
    snap_files = sim["snap_files"]
    snap_numbers = sim["numbers"]

    time_data = timeSeriesData(snap_files; sim_cosmo)
    time_unit = time_data["units"]["time"]
    clock = time_data["clock_time"]

    # Create a directory to store the figures, if it doesn't exist.
    mkpath(output_path * "images/")

    # Generate and store the plots.
    short_snaps = @view snap_files[1:step:end] 
    
    # animation = @animate          
    for (i, snapshot) in enumerate(short_snaps)

        temp_data = tempData(snapshot; sim_cosmo, temp_unit)
        density_data = densityData(snapshot; sim_cosmo, density_unit)

        # Snashot number.
        number = snap_numbers[1 + step * (i - 1)]

        savefig(
            rhoTempPlot(
                temp_data,
                density_data::Dict{String, Any},
                clock[1 + step * (i - 1)] * time_unit,
            ),
            output_path * "images/" * base_name * "_" * number * format,
        )
        
    end

    # Make the GIF.
    # gif(animation, output_path * anim_name * ".gif", fps = frame_rate)

    # Make the video.
    # makeVideo(output_path * "images/", format, output_path, anim_name, frame_rate)

    return nothing
end

"""
    KennicuttSchmidtPipeline(
        base_name::String,
        source_path::String; 
        <keyword arguments>
    )::Nothing

Save the results of the KennicuttSchmidtPlot function as one image per snapshot.

It will produce output only for the snapshots that have enough young stars to produce 
at least five data points for the fitting.

# Arguments
- `base_name::String`: Base names of the snapshot files, set in the GADGET 
  variable SnapshotFileBase.
- `source_path::String`: Paths to the directories containing the snapshot files, 
  set in the GADGET variable OutputDir.
- `output_path::String = "Kennicutt_Schmidt/"`: Path to the output directory. The images 
  will be stored in `output_path`images/ and will be named `base_name`_XXX`format` where XXX 
  is the number of the snapshot. The GIF and the video will be stored in `output_path`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `step::Int64 = 1`: Step used to traverse the list of snapshots. The default is 1, 
  i.e. all snapshots will be plotted.
- `temp_filter::Unitful.Quantity`: Maximum temperature allowed for the gas particles.
- `age_filter::Unitful.Quantity`: Maximum age allowed for the stars.
- `box_size::Unitful.Quantity = 1000UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. It has to have units, e.g. 1000UnitfulAstro.kpc, 
  which is the default. Its units don't have to be the same as `length_unit`.
- `bins::Int64 = 50`: Number of subdivisions of [0, `max_r`] to be used. 
  It has to be at least 5.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Myr, which is the default.
- `mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun`: Unit of mass to be used in the output, 
  all available mass units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.Msun, which is the default.
- `temp_unit::Unitful.FreeUnits = Unitful.K`: Unit of temperature to be used in the 
  output, all available temperature units in Unitful and UnitfulAstro can be used, 
  e.g. Unitful.K, which is the default.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used in the 
  output, all available length units in Unitful and UnitfulAstro can be used, 
  e.g. UnitfulAstro.kpc, which is the default.
- `format::String = ".png"`: File format of the output figure. All formats supported by 
  pgfplotsx can be used, namely ".pdf", ".tex", ".svg" and ".png", which is the default. 
"""
function KennicuttSchmidtPipeline(
    base_name::String,
    source_path::String;
    output_path::String = "Kennicutt_Schmidt/",
    sim_cosmo::Int64 = 0,
    step::Int64 = 1,
    temp_filter::Unitful.Quantity = 3e4Unitful.K,
    age_filter::Unitful.Quantity = 200UnitfulAstro.Myr,
    box_size::Unitful.Quantity = 1000UnitfulAstro.kpc,
    bins::Int64 = 50,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun,
    temp_unit::Unitful.FreeUnits = Unitful.K,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    format::String = ".png",
)::Nothing

    # Get the simulation data.
    sim = getSnapshots(base_name, source_path)
    snap_files = sim["snap_files"]
    snap_numbers = sim["numbers"]

    time_data = timeSeriesData(snap_files; sim_cosmo)
    clock_unit = time_data["units"]["time"]
    # Clock time of each snapshot.
    clock = time_data["clock_time"]

    # Create a directory to store the figures, if it doesn't exist.
    mkpath(output_path * "images/")

    # Generate and store the plots.
    short_snaps = @view snap_files[1:step:end] 

    @inbounds for (i, snapshot) in enumerate(short_snaps)

        header = read_header(snapshot)
        if header.nall[5] != 0
            
            number = snap_numbers[1 + step * (i - 1)]
            # Clock time of the snapshot i.
            now = uconvert(time_unit, clock[1 + step * (i - 1)] * clock_unit)

            # Gas masses.
            gas_mass_data = massData(snapshot, "gas"; sim_cosmo, mass_unit)
            # Gas temperatures.
            temperature_data = tempData(snapshot; sim_cosmo, temp_unit)
            # Stars masses.
            star_mass_data = massData(snapshot, "stars"; sim_cosmo, mass_unit)
            # Stars ages.
            age_data = ageData(snapshot, now; sim_cosmo)
            # Positions.
            pos_data = positionData(snapshot; sim_cosmo, box_size, length_unit)

            figure = KennicuttSchmidtPlot(
                gas_mass_data,
                temperature_data,
                star_mass_data,
                age_data,
                pos_data,
                temp_filter,
                age_filter,
                box_size,
                clock[1 + step * (i - 1)] * clock_unit;
                bins,
            )

            if figure !== nothing
                # If there was enough data to make a fit.
                savefig(
                    figure,
                    output_path * "images/" * base_name * "_" * number * format,
                )
            end

        end
    end

    return nothing
end