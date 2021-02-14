############################################################################################
# DATA ACQUISITION FUNCTIONS.
############################################################################################

"""
    function getSnapshots(
        base_name::String,
        source_path::String,
    )::Dict{String, Vector{String}}

Give the paths to the GADGET output files, grouping them by snapshot in a Tuple.

# Arguments
- `base_name::String`: Base name of the snapshot files, 
  set in the GADGET variable SnapshotFileBase.
- `source_path::String`: Path to the directory containing the snapshot files, 
  set in the GADGET variable OutputDir.

# Returns
- A dictionary with two entries. 
  - Key "numbers" => A tuple of Strings with the numbers 
  that characterize each snapshot.
  - Key "snap_files" => A tuple of Strings with the paths to the snapshot files,
  named as follows: (`source_path` `base_name`_001 ... `source_path` `base_name`_NNN) 
  where NNN is the total number of snapshots.
"""
function getSnapshots(base_name::String, source_path::String)::Dict{String, Vector{String}}

    # Get the full list of paths to every GADGET file in `source_path`.
    file_list = glob(base_name * "_*", source_path)

    # Data availability check.
    !isempty(file_list) || error("I couldn't find any snapshots in $source_path.")

    # Get the number of files per snapshot.
    num_files = read_header(first(file_list)).num_files

    if num_files > 1
        # If there are multiple files per snapshot, delete the trailing '.n'.
        map!(x -> rsplit(x, "."; limit = 2)[1], file_list, file_list)
        # And delete duplicates.
        unique!(file_list)
    end

    # Get the numbers that characterize each snapshot.
    numbers = map(x -> rsplit(x, base_name * '_'; limit = 2)[2], file_list)

    return Dict("numbers" => numbers, "snap_files" => file_list)
end

"""
    timeSeriesData(snap_files::Vector{String}; <keyword arguments>)::Dict{String, Any}
					
Get several parameters defined at every snapshot, as a series of values for the whole 
simulation. 

The parameters are:

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
- "gas_density" (Global gas density)	               
- "gas_frac" (Gas fraction)		                
- "dm_frac" (Dark matter fraction)		                
- "star_frac" (Star fraction)	                   
- "gas_bar_frac" (Baryonic gas fraction)                  
- "star_bar_frac" (Baryonic star fraction)

# Arguments
- `snap_files::Vector{String}`: Output of the function getSnapshots corresponding 
  to the key "snap_files", containing an Array with the paths to the snapshots.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun`: Unit of mass to be used in the output, 
  all available mass units in Unitful.jl and UnitfulAstro.jl can be used, 
  e.g. UnitfulAstro.Msun, which is the default.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in Unitful.jl and UnitfulAstro.jl can be used, 
  e.g. UnitfulAstro.Myr, which is the default.
- `sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr`: Unit of mass/time to 
  be used in the output, all available time and mass units in Unitful.jl and UnitfulAstro.jl 
  can be used, e.g. UnitfulAstro.Msun / UnitfulAstro.yr, which is the default.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used 
  in the output, all available length units in Unitful.jl and UnitfulAstro.jl 
  can be used, e.g. UnitfulAstro.kpc, which is the default.

# Returns
- A dictionary.
  - Key "{property}" => One dimensional array with the numeric values 
  of the {property} in the key (one value per snapshot) for the whole simulation. 
  - Key "units" => A dictionary with the units used, for easy piping with other functions.
  - Key "labels" => A dictionary with the labels to be used when plotting the quantities.
"""
function timeSeriesData(
    snap_files::Vector{String};
    sim_cosmo::Int64 = 0,
    mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
)::Dict{String, Any}

    # Number of sanpshots.
    n_files = length(snap_files)

    # Output data structure.
    time_series = Dict(
        "scale_factor" => Vector{Float64}(undef, n_files),  # Dimensionless.
        "redshift" => Vector{Float64}(undef, n_files),      # Dimensionless.
        "clock_time" => Vector{Float64}(undef, n_files),    # `time_unit`.
        "sfr" => Vector{Float64}(undef, n_files),           # `sfr_unit`.
        "sfr_prob" => Vector{Float64}(undef, n_files),      # `sfr_unit`.
        "gas_number" => Vector{Int64}(undef, n_files),      # Dimensionless.
        "dm_number" => Vector{Int64}(undef, n_files),       # Dimensionless.
        "star_number" => Vector{Int64}(undef, n_files),     # Dimensionless.
        "gas_mass" => Vector{Float64}(undef, n_files),      # `mass_unit`.
        "dm_mass" => Vector{Float64}(undef, n_files),       # `mass_unit`.
        "star_mass" => Vector{Float64}(undef, n_files),     # `mass_unit`.
        "gas_density" => Vector{Float64}(undef, n_files),   # `mass_unit` / `length_unit`^3.
        "gas_frac" => Vector{Float64}(undef, n_files),      # Dimensionless.
        "dm_frac" => Vector{Float64}(undef, n_files),       # Dimensionless.
        "star_frac" => Vector{Float64}(undef, n_files),     # Dimensionless.
        "gas_bar_frac" => Vector{Float64}(undef, n_files),  # Dimensionless.
        "star_bar_frac" => Vector{Float64}(undef, n_files), # Dimensionless.
        # Unit pass-through.
        "units" => Dict(
            "mass" => mass_unit,
            "time" => time_unit,
            "sfr" => sfr_unit,
            "length" => length_unit,
        ),
        # Labels for printing.
        "labels" => Dict(
            "scale_factor" => "a",
            "redshift" => "z",
            "clock_time" => "t",
            "sfr" => "SFR",
            "sfr_prob" => "SFR probability",
            "gas_number" => "Gas particle number",
            "dm_number" => "Dark matter particle number",
            "star_number" => "Star number",
            "gas_mass" => "Total gas mass",
            "dm_mass" => "Total dark matter mass",
            "star_mass" => "Total star mass",
            "gas_density" => "Total gas density",
            "gas_frac" => "Gas fraction",
            "dm_frac" => "Dark matter fraction",
            "star_frac" => "Star fraction",
            "gas_bar_frac" => "Baryonic gas fraction",
            "star_bar_frac" => "Baryonic star fraction",
        ),
    )

    for (i, snapshot) in enumerate(snap_files)

        # Data availability check.
        (
            block_present(snapshot, "MASS") ||
            error("There is no block 'MASS' in snapshot located at $snapshot")
        )
        (
            block_present(snapshot, "SFR") ||
            error("There is no block 'SFR' in snapshot located at $snapshot")
        )

        header = read_header(snapshot)

        # Struct for unit conversion.
        if sim_cosmo == 1
            GU = GadgetPhysicalUnits(a_scale = header.time, hpar = header.h0)

            a = header.time     # Scale factor.
            z = (1 / a) - 1     # Redshift.

            # Ω_K (curvature)
            omega_K = 1 - header.omega_0 - header.omega_l
            # Energy function.
            E = header.omega_0 / (a * a * a) + omega_K / (a * a) + header.omega_l
            # Hubble constant in 1 / Gyr.
            H = header.h0 * HUBBLE_CONST * sqrt(E)
            # Integrand of the time integral in Gyr. 
            t_int = 1 / (H * a)
        else
            # For Newtonian simulations uses the default scale factor: a = 1.
            GU = GadgetPhysicalUnits(hpar = header.h0)

            a = 1.0      # Scale factor.
            z = 0.0     # Redshift.
        end

        # Number of particles of each type.
        gas_number = header.nall[1]
        dm_number = header.nall[2]
        star_number = header.nall[5]

        # Total mass of gas in `mass_unit`.
        if gas_number != 0
            if header.massarr[1] != 0
                # If all gas particles have the same mass.
                masses = fill(header.massarr[1], (gas_number))
                gas_mass = masses[1] * gas_number
            else
                masses = read_snap(snapshot, "MASS", 0)
                gas_mass = sum(masses)
            end
            gas_mass = ustrip(Float64, mass_unit, gas_mass * GU.m_cgs)

            # Global gas density.
            densities = densityData(
                snapshot; 
                sim_cosmo, 
                density_unit = mass_unit / length_unit^3
            )
            volume = sum(masses ./ densities["gas_density"])
            gas_density = gas_mass / volume
        else
            # In the case that there are no gas particles.
            gas_mass = 0.0
            gas_density = 0.0
        end

        # Total mass of dark matter in `mass_unit`.
        if dm_number != 0
            if header.massarr[2] != 0
                # If all dark matter particles have the same mass.
                dm_mass = header.massarr[2] * dm_number
            else
                dm_mass = sum(read_snap(snapshot, "MASS", 1))
            end
            dm_mass = ustrip(Float64, mass_unit, dm_mass * GU.m_cgs)
        else
            # In the case that there are no dark matter particles.
            dm_mass = 0.0
        end

        # Total mass of stars in `mass_unit`.
        if star_number != 0
            if header.massarr[5] != 0
                # If all star particles have the same mass.
                star_mass = header.massarr[5] * star_number
            else
                star_mass = sum(read_snap(snapshot, "MASS", 4))
            end
            star_mass = ustrip(Float64, mass_unit, star_mass * GU.m_cgs)
        else
            # In the case that there are no star particles.
            star_mass = 0.0
        end

        # Total baryonic mass in `mass_unit`.
        baryonic_mass = gas_mass + star_mass
        # Total mass of the system in `mass_unit`.
        total_mass = baryonic_mass + star_mass

        # Physical SFR and time in `sfr_unit` and `time_unit` respectively.
        if i == 1
            # In the first step, clock time and SFR are set to 0.
            t = 0.0
            sfr = 0.0
        else
            # In every other step, clock time and SFR are calculated. 
            if sim_cosmo == 1
                # Integrand of the time integral in `time_unit`.
                t_int = ustrip(Float64, time_unit, t_int * UnitfulAstro.Gyr)
                Δa = a - time_series["scale_factor"][i - 1]
                Δt = Δa * t_int
                Δstar_mass = star_mass - time_series["star_mass"][i - 1]

                t = time_series["clock_time"][i - 1] + Δt
                sfr = ustrip(Float64, sfr_unit, (Δstar_mass / Δt) * (mass_unit / time_unit))
            else
                t = ustrip(Float64, time_unit, header.time * GU.t_Myr)

                Δt = t - time_series["clock_time"][i - 1]
                Δstar_mass = star_mass - time_series["star_mass"][i - 1]

                sfr = ustrip(Float64, sfr_unit, (Δstar_mass / Δt) * (mass_unit / time_unit))
            end
        end

        # SFR given by the classic prescription of GADGET in `sfr_unit`.
        if header.nall[1] != 0
            # Conversion factor from internal units of time to physical units of time.
            if sim_cosmo == 1
                t_conv = t_int * UnitfulAstro.Gyr
            else
                t_conv = GU.t_Myr
            end

            sfr_prob = sum(read_snap(snapshot, "SFR", 0))
            sfr_prob = ustrip(Float64, sfr_unit, sfr_prob * (GU.m_cgs / t_conv))
        else
            # In the case that there are no gas particles.
            sfr_prob = 0.0
        end

        # Time-like parameters.
        time_series["scale_factor"][i] = a
        time_series["redshift"][i] = z
        time_series["clock_time"][i] = t

        # Real SFR.
        time_series["sfr"][i] = sfr
        # SFR given by the classic prescription (∝ probability of star formation).
        time_series["sfr_prob"][i] = sfr_prob

        # Number of particles of each type.
        time_series["gas_number"][i] = gas_number
        time_series["dm_number"][i] = dm_number
        time_series["star_number"][i] = star_number

        # Total mass of each type of particle.
        time_series["gas_mass"][i] = gas_mass
        time_series["dm_mass"][i] = dm_mass
        time_series["star_mass"][i] = star_mass

        # Global gas density.
        time_series["gas_density"][i] = gas_density

        # Mass fraction relative to the total mass of the system.
        time_series["gas_frac"][i] = gas_mass / total_mass
        time_series["dm_frac"][i] = dm_mass / total_mass
        time_series["star_frac"][i] = star_mass / total_mass

        # Mass fraction relative to the total baryonic mass.
        time_series["gas_bar_frac"][i] = gas_mass / baryonic_mass
        time_series["star_bar_frac"][i] = star_mass / baryonic_mass
    end

    return time_series
end

"""
    positionData(snapshot::String; <keyword arguments>)::Dict{String, Any}

Get the coordinates of the particles at a specific time step.

# Arguments
- `snapshot::String`: Path to a given snapshot.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `box_size::Unitful.Quantity = 1000UnitfulAstro.kpc`: Size of the plotting region 
  if vacuum boundary conditions were used. It has to have units, e.g. 1000UnitfulAstro.kpc, 
  which is the default. Its units don't have to be the same as `length_unit`.
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used 
  in the output, all available length units in Unitful.jl and UnitfulAstro.jl 
  can be used, e.g. UnitfulAstro.kpc, which is the default.

# Returns
- A dictionary with six entries.
  - Keys "gas", "dark_matter", "stars" => 2 dimensional arrays which contains 
  the positions of the particles of the type given by the key. Each row is a 
  particle and each column corresponds to coordinates x, y and z respectively.
  - Key "box_size" => The range of values for the plotting of the positions, 
  i.e. a range of ± `box_size` if vacuum boundary conditions were used, 
  or (0, `header.boxsize`) if periodic boundary conditions were used.
  Notice how the side length of the region is 2 * `box_size` for vacuum boundary 
  conditions and `header.boxsize` for periodic boundary conditions.
  - Key "periodic" => A Boolean indicating the type of boundary condition.
  false -> vacuum boundary condition.
  true -> periodic boundary condition.
  - Key "unit" => The unit of lengthused, i.e. is a pass-through of `length_unit`. 
"""
function positionData(
    snapshot::String;
    sim_cosmo::Int64 = 0,
    box_size::Unitful.Quantity = 1000UnitfulAstro.kpc,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
)::Dict{String, Any}

    # Data availability check.
    (
        block_present(snapshot, "POS") ||
        error("There is no block 'POS' in snapshot located at $snapshot")
    )

    header = read_header(snapshot)

    # Struct for unit conversion.
    if sim_cosmo == 1
        GU = GadgetPhysicalUnits(a_scale = header.time, hpar = header.h0)
    else
        # For Newtonian simulation uses the default scale factor: a = 1.
        GU = GadgetPhysicalUnits(hpar = header.h0)
    end

    # Get the correct region size in `length_unit`, 
    # given the type of boundary condition used.
    if header.boxsize == 0
        # Vacuum boundary condition.
        size = round(ustrip(Float64, length_unit, box_size), sigdigits = 1)
    else
        # Periodic boundary conditions.
        size = round(ustrip(Float64, length_unit, header.boxsize * GU.x_cgs), sigdigits = 1)
    end

    if header.nall[1] != 0
        gas_pos = read_snap(snapshot, "POS", 0)
        # Transformation from internal units to `length_unit`.
        gas_pos = @. ustrip(Float64, length_unit, gas_pos * GU.x_cgs)
    else
        # In the case that there are no gas particles.
        gas_pos = Array{Float64}(undef, 3, 0)
    end

    if header.nall[2] != 0
        dm_pos = read_snap(snapshot, "POS", 1)
        # Transformation from internal units to `length_unit`.
        dm_pos = @. ustrip(Float64, length_unit, dm_pos * GU.x_cgs)
    else
        # In the case that there are no dark matter particles.
        dm_pos = Array{Float64}(undef, 3, 0)
    end

    if header.nall[5] != 0
        star_pos = read_snap(snapshot, "POS", 4)
        # Transformation from internal units to `length_unit`.
        star_pos = @. ustrip(Float64, length_unit, star_pos * GU.x_cgs)
    else
        # In the case that there are no star particles.
        star_pos = Array{Float64}(undef, 3, 0)
    end

    return Dict(
        "gas" => gas_pos,
        "dark_matter" => dm_pos,
        "stars" => star_pos,
        "box_size" => size,
        "periodic" => (header.boxsize != 0),
        "unit" => length_unit,
    )
end

"""
    densityData(snapshot::String; <keyword arguments>)::Dict{String,Any}

Get the densities of the gas particles at a specific time step.

# Arguments
- `snapshot::String`: Path to a given snapshot.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `density_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.kpc^3`: Unit of density 
  to be used in the output, all available density units in Unitful.jl and UnitfulAstro.jl can 
  be used, e.g. UnitfulAstro.Msun / UnitfulAstro.kpc^3, which is the default.

# Returns
- A dictionary with two entries.
  - Key "gas_density" => One dimensional array, which contains the densities
    of the gas particles. 
  - Key "unit" => The unit of density used, i.e. is a pass-through of `density_unit`. 
"""
function densityData(
    snapshot::String;
    sim_cosmo::Int64 = 0,
    density_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.kpc^3,
)::Dict{String, Any}

    # Data availability check.
    (
        block_present(snapshot, "RHO") ||
        error("There is no block 'RHO' in snapshot located at $snapshot")
    )

    header = read_header(snapshot)

    # Struct for unit conversion.
    if sim_cosmo == 1
        GU = GadgetPhysicalUnits(a_scale = header.time, hpar = header.h0)
    else
        # For Newtonian simulation uses the default scale factor: a = 1.
        GU = GadgetPhysicalUnits(hpar = header.h0)
    end

    if header.nall[1] != 0
        ρ = read_snap(snapshot, "RHO", 0)
        # Transformation from internal units to `density_unit`.
        ρ = @. ustrip(Float64, density_unit, ρ * GU.rho_cgs)
    else
        # In the case that there are no gas particles.
        ρ = Float64[]
    end

    return Dict("gas_density" => ρ, "unit" => density_unit)
end

"""
    hsmlData(snapshot::String; <keyword arguments>)::Dict{String,Any}

Get the smoothing lengths of the gas particles at a specific time step.

# Arguments
- `snapshot::String`: Path to a given snapshot.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `length_unit::Unitful.FreeUnits = UnitfulAstro.kpc`: Unit of length to be used 
  in the output, all available length units in Unitful.jl and UnitfulAstro.jl 
  can be used, e.g. UnitfulAstro.kpc, which is the default.

# Returns
- A dictionary with two entries.
    - Key "gas_hsml" => A one dimensional array, which contains the smoothing lengths
    of the particles of gas. 
    - Key "unit" => The unit used, i.e. is a pass-through of `length_unit`. 
"""
function hsmlData(
    snapshot::String;
    sim_cosmo::Int64 = 0,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
)::Dict{String, Any}

    # Data availability check.
    (
        block_present(snapshot, "HSML") ||
        error("There is no block 'HSML' in snapshot located at $snapshot")
    )

    header = read_header(snapshot)

    # Struct for unit conversion.
    if sim_cosmo == 1
        GU = GadgetPhysicalUnits(a_scale = header.time, hpar = header.h0)
    else
        # For Newtonian simulation uses the default scale factor: a = 1.
        GU = GadgetPhysicalUnits(hpar = header.h0)
    end

    if header.nall[1] != 0
        hsml = read_snap(snapshot, "HSML", 0)
        # Transformation from internal units to `length_unit`.
        hsml = @. ustrip(Float64, length_unit, hsml * GU.x_cgs)
    else
        # In the case that there are no gas particles.
        hsml = Float64[]
    end

    return Dict("gas_hsml" => hsml, "unit" => length_unit)
end

"""
    massData(snapshot::String, type::String; <keyword arguments>)::Dict{String,Any}

Get the mass of the particles at a specific time step.

# Arguments
- `snapshot::String`: Path to a given snapshot.
- `type::String`: Particle type.
  "gas" -> Gas particle. 
  "dark_matter" -> Dark matter particle.
  "stars" -> Star particle.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun`: Unit of mass to be used in the output, 
  all available mass units in Unitful.jl and UnitfulAstro.jl can be used, 
  e.g. UnitfulAstro.Msun, which is the default.

# Returns
- A dictionary with three entries.
  - Key "mass" => A one dimensional array which contains the masses of the particles. 
  - Key "unit" => The unit of mass used, i.e. is a pass-through of `mass_unit`. 
  - Key "type" => A String with the particle type.
"""
function massData(
    snapshot::String,
    type::String;
    sim_cosmo::Int64 = 0,
    mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun,
)::Dict{String, Any}

    # Data availability check.
    (
        block_present(snapshot, "MASS") ||
        error("There is no block 'MASS' in snapshot located at $snapshot")
    )

    header = read_header(snapshot)

    # Struct for unit conversion.
    if sim_cosmo == 1
        GU = GadgetPhysicalUnits(a_scale = header.time, hpar = header.h0)
    else
        # For Newtonian simulation uses the default scale factor: a = 1.
        GU = GadgetPhysicalUnits(hpar = header.h0)
    end

    # Select type of particle.
    if type == "gas"
        type_num = 0
    elseif type == "dark_matter"
        type_num = 1
    elseif type == "stars"
        type_num = 4
    else
        error("Particle type '$type' not supported. 
        The supported types are 'gas', 'dark_matter' and 'stars'")
    end

    if header.nall[type_num + 1] != 0
        m = read_snap(snapshot, "MASS", type_num)
        # Transformation from internal units to `mass_unit`.
        m = @. ustrip(Float64, mass_unit, m * GU.m_cgs)
    else
        # In the case that there are no particles.
        m = [Inf]
    end

    return Dict("mass" => m, "unit" => mass_unit, "type" => type)
end

"""
    zData(snapshot::String, type::String; <keyword arguments>)::Dict{String,Any}

Get the metallicity (as mass content of metals) of the particles at a specific time step.

# Arguments
- `snapshot::String`: Path to a given snapshot.
- `type::String`: Particle type.
  "gas" -> Gas particle. 
  "stars" -> Star particle.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun`: Unit of mass to be used in the output, 
  all available mass units in Unitful.jl and UnitfulAstro.jl can be used, 
  e.g. UnitfulAstro.Msun, which is the default.

# Returns
- A dictionary with two entries.
  - Key "Z" => A one dimensional array which contains the metallicities of the particles.  
  - Key "unit" => The unit of mass used, i.e. is a pass-through of `mass_unit`. 
  - Key "type" => A String with the particle type. 
"""
function zData(
    snapshot::String,
    type::String;
    sim_cosmo::Int64 = 0,
    mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun,
)::Dict{String, Any}

    # Data availability check.
    (
        block_present(snapshot, "Z") ||
        error("There is no block 'Z' in snapshot located at $snapshot")
    )

    header = read_header(snapshot)

    # Struct for unit conversion.
    if sim_cosmo == 1
        GU = GadgetPhysicalUnits(a_scale = header.time, hpar = header.h0)
    else
        # For Newtonian simulation uses the default scale factor: a = 1.
        GU = GadgetPhysicalUnits(hpar = header.h0)
    end

    # Select type of particle.
    if type == "gas"
        type_num = 0
    elseif type == "stars"
        type_num = 4
    else
        error("Particle type '$type' not supported. 
        The supported types are 'gas' and 'stars'")
    end

    if header.nall[type_num + 1] != 0
        z = read_snap(snapshot, "Z", type_num)

        # Initialize output array.
        Z = similar(Array{Float64}, axes(z, 2))
        for i in eachindex(Z)
            # Add up all elements but the ones at position 1 and 7, i.e. H and He.
            z_tot = sum(z[[2, 3, 4, 5, 6, 8, 9, 10, 11, 12], i])
            # Transformation from internal units to `mass_unit`.
            z_tot = ustrip(Float64, mass_unit, z_tot * GU.m_cgs)

            @inbounds Z[i] = z_tot
        end
    else
        # In the case that there are no particles.
        Z = [Inf]
    end

    return Dict("Z" => Z, "unit" => mass_unit, "type" => type)
end

"""
    tempData(snapshot::String; <keyword arguments>)::Dict{String,Any}

Get the temperature of the gas particles at a specific time step.

# Arguments
- `snapshot::String`: Path to a given snapshot.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `temp_unit::Unitful.FreeUnits = Unitful.K`: Unit of temperature to be used in the output, 
  all available temperature units in Unitful.jl and UnitfulAstro.jl can be used, 
  e.g. Unitful.K, which is the default.

# Returns
- A dictionary with two entries.
  - Key "T" => A one dimensional array which contains the temperatures of the particles.  
  - Key "unit" => The unit of temperature used, i.e. is a pass-through of `temp_unit`. 
"""
function tempData(
    snapshot::String;
    sim_cosmo::Int64 = 0,
    temp_unit::Unitful.FreeUnits = Unitful.K,
)::Dict{String, Any}

    # Data availability check.
    (
        block_present(snapshot, "U") ||
        error("There is no block 'U' in snapshot located at $snapshot")
    )
    (
        block_present(snapshot, "NE") ||
        error("There is no block 'NE' in snapshot located at $snapshot")
    )
    (
        block_present(snapshot, "Z") ||
        error("There is no block 'Z' in snapshot located at $snapshot")
    )

    header = read_header(snapshot)

    # Struct for unit conversion.
    if sim_cosmo == 1
        GU = GadgetPhysicalUnits(a_scale = header.time, hpar = header.h0)
    else
        # For Newtonian simulation uses the default scale factor: a = 1.
        GU = GadgetPhysicalUnits(hpar = header.h0)
    end

    if header.nall[1] != 0

        # Mass.
        mass_data = massData(snapshot, "gas"; sim_cosmo)
        mass = mass_data["mass"]
        mass_unit = mass_data["unit"]

        # Metallicity.
        z = read_snap(snapshot, "Z", 0)
        Z = @. ustrip(Float64, mass_unit, z * GU.m_cgs)

        # Internal enery per unit mass.
        u = read_snap(snapshot, "U", 0)
        U = @. uconvert(Unitful.J / mass_unit, u * (GU.E_cgs / GU.m_cgs))

        # ne := number_of_electrons / number_of_Hydrogen_atoms.
        ne = read_snap(snapshot, "NE", 0)

        # xH := mass_fraction_of_Hydrogen.
        xH = Z[7, :] ./ mass 
        # yHe := number_of_Helium_atoms / number_of_Hydrogen_atoms.
        yHe = @. (1.0 - xH) / (4.0 * xH)
        # μ := total_mass / (total_number_of_particles * proton_mass).
        μ = @. (1.0 + 4.0 * yHe) / (1.0 + yHe + ne)
        # T = (adiabatic_index - 1) * internal_energy_per_unit_mass * 
        #     (total_mass / total_number_of_particles) / Boltzmann_constant
        T = @. ustrip(Float64, temp_unit, 2/3 * U * μ * Unitful.mp / Unitful.k)
        
    else
        # In the case that there are no particles.
        T = [Inf]
    end

    return Dict("T" => T, "unit" => temp_unit)
end

"""
    ageData(snapshot::String, now::Float64; <keyword arguments>)::Dict{String,Any}

Get the ages of the stars at a specific time step.

# Arguments
- `snapshot::String`: Path to a given snapshot.
- `now::Unitful.Quantity`: Clock time of `snapshot`, with units. All available time units in 
  Unitful.jl and UnitfulAstro.jl can be used.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).

# Returns
- A dictionary with two entries.
  - Key "ages" => The ages of the stars.  
  - Key "unit" => The unit of time used. 
"""
function ageData(
    snapshot::String,
    now::Unitful.Quantity;
    sim_cosmo::Int64 = 0,
)::Dict{String, Any}

    # Data availability check.
    (
        block_present(snapshot, "AGE") ||
        error("There is no block 'AGE' in snapshot located at $snapshot")
    )

    header = read_header(snapshot)

    # Struct for unit conversion.
    if sim_cosmo == 1
        GU = GadgetPhysicalUnits(a_scale = header.time, hpar = header.h0)

        a = header.time     # Scale factor.
        z = (1 / a) - 1     # Redshift.

        # Ω_K (curvature)
        omega_K = 1 - header.omega_0 - header.omega_l
        # Energy function.
        E = header.omega_0 / (a * a * a) + omega_K / (a * a) + header.omega_l
        # Hubble constant in 1 / Gyr.
        H = header.h0 * HUBBLE_CONST * sqrt(E)
        # Integrand of the time integral in Gyr. 
        t_int = 1 / (H * a)
        # Conversion factor from internal units of time to physical units of time.
        t_conv = t_int * UnitfulAstro.Gyr
    else
        # For Newtonian simulation uses the default scale factor: a = 1.
        GU = GadgetPhysicalUnits(hpar = header.h0)

        # Conversion factor from internal units of time to physical units of time.
        t_conv = GU.t_Myr
    end

    # Stars time of birth. 
    birth_time = read_snap(snapshot, "AGE", 4)

    # Unit conversion.
    time_unit = unit(now)
    birth_time = @. ustrip(Float64, time_unit, birth_time * t_conv)
    
    # Stars ages. 
    clock_now = ustrip(now)
    ages = clock_now .- birth_time

    return Dict("ages" => ages, "unit" => time_unit)
end

"""
    birthPlace(
        snap_index::Int64,
        snap_files::Vector{String},
        time_stamps::Vector{Float64},
        stamps_unit::Unitful.FreeUnits; 
        <keyword arguments>
    )::Array{Float64, 2}

Get the birth location of the stars in a given snapshot.

# Arguments
- `snap_index::Int64`: Index in `snap_files` of the snapshot whose stars will be located.
- `snap_files::Vector{String}`: Output of the function getSnapshots corresponding 
  to the key "snap_files", containing an Array with the paths to the snapshots.
- `time_stamps::Vector{Float64}`: Clock time of every snapshot in `snap_files`.
- `stamps_unit::Unitful.FreeUnits`: Unit of time of the `time_stamps`.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time of `time_stamps`.

# Returns
-  2 dimensional arrays which contains the positions of the stars. Each row is a star
  and each column corresponds to coordinates x, y and z respectively.
"""
function birthPlace(
    snap_index::Int64,
    snap_files::Vector{String},
    time_stamps::Vector{Float64},
    stamps_unit::Unitful.FreeUnits;
    sim_cosmo::Int64 = 0,
    length_unit::Unitful.FreeUnits = UnitfulAstro.kpc,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
)::Dict{String, Any}

    # Bounds check.
    (
        1 < snap_index <= length(snap_files) ||
        throw(BoundsError("There is no snapshot and index $snap_index."))
    )
    # Dimension consistency check.
    (
        length(time_stamps) == length(snap_files) ||
        throw(DimensionMismatch("The input vectors should have the same length."))
    )

    snapshot = snap_files[snap_index]

    # Data availability check.
    (
        block_present(snapshot, "ID") ||
        error("There is no block 'ID' in snapshot located at $snapshot")
    )
    (
        block_present(snapshot, "AGE") ||
        error("There is no block 'AGE' in snapshot located at $snapshot")
    )

    header = read_header(snapshot)

    # Struct for unit conversion.
    if sim_cosmo == 1
        GU = GadgetPhysicalUnits(a_scale = header.time, hpar = header.h0)

        a = header.time     # Scale factor.
        z = (1 / a) - 1     # Redshift.

        # Ω_K (curvature)
        omega_K = 1 - header.omega_0 - header.omega_l
        # Energy function.
        E = header.omega_0 / (a * a * a) + omega_K / (a * a) + header.omega_l
        # Hubble constant in 1 / Gyr.
        H = header.h0 * HUBBLE_CONST * sqrt(E)
        # Integrand of the time integral in Gyr. 
        t_int = 1 / (H * a)
        # Conversion factor from internal units of time to physical units of time.
        t_conv = t_int * UnitfulAstro.Gyr
    else
        # For Newtonian simulation uses the default scale factor: a = 1.
        GU = GadgetPhysicalUnits(hpar = header.h0)

        # Conversion factor from internal units of time to physical units of time.
        t_conv = GU.t_Myr
    end

    # Stars time of birth. 
    ages = read_snap(snapshot, "AGE", 4)
    birth_times = @. ustrip(Float64, time_unit, ages * t_conv)
    time_stamps = @. ustrip(Float64, time_unit, time_stamps * stamps_unit)
    # Stars IDs.
    ids = read_snap(snapshot, "ID", 4)

    birth_place = Vector{Float64}[]
    for (id, birth_time) in zip(ids, birth_times)

        # Index of the snapshot where the target star was born.
        snap_idx = findfirst(x -> x >= birth_time, time_stamps)

        # IDs of the stars in the snapshot where the target star was born.
        nursery_ids = read_snap(snap_files[snap_idx], "ID", 4)

        # Index of the target star, in the snapshot where it was born.
        birth_idx = findfirst(x -> x == id, nursery_ids)
        while birth_idx === nothing
            snap_idx += 1
            nursery_ids = read_snap(snap_files[snap_idx], "ID", 4)
            birth_idx = findfirst(x -> x == id, nursery_ids)
        end

        # Position of the target star, in the snapshot where it was born.
        raw_pos = read_snap(snap_files[snap_idx], "POS", 4)[:, birth_idx]
        nursery_pos = @. ustrip(Float64, length_unit, raw_pos * GU.x_cgs)

        push!(birth_place, nursery_pos)
    end

    return Dict("birth_place" => hcat(birth_place...), "unit" => length_unit)
end

"""
    sfrTxtData(
        source_path::String,
        base_name::String; 
        <keyword arguments>
    )::Dict{Int64, Vector{Float64}}

Get the column data from the file sfr.txt.

Transform from internal units to the ones given by `mass_unit`, `time_unit` and `sfr_unit`.

# Warning 
THIS FUNCTION NEEDS A FILE (sfr.txt) WHICH IS NOT PART OF ANY PUBLIC VERSION OF GADGET.

# Arguments
- `source_path::String`: Path to the directory containing the sfr.txt file.
- `base_name::String`: Base name of the snapshot files, 
  set in the GADGET variable SnapshotFileBase.
- `sim_cosmo::Int64 = 0`: Value of the GADGET variable ComovingIntegrationOn: 
  0 -> Newtonian simulation (static universe).
  1 -> Cosmological simulation (expanding universe).
- `mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun`: Unit of mass to be used in the output, 
  all available mass units in Unitful.jl and UnitfulAstro.jl can be used, 
  e.g. UnitfulAstro.Msun, which is the default.
- `time_unit::Unitful.FreeUnits = UnitfulAstro.Myr`: Unit of time to be used in the output, 
  all available time units in Unitful.jl and UnitfulAstro.jl can be used, 
  e.g. UnitfulAstro.Myr, which is the default.
- `sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr`: Unit of mass/time to 
  be used in the output, all available time and mass units in Unitful.jl and UnitfulAstro.jl 
  can be used, e.g. UnitfulAstro.Msun/UnitfulAstro.yr, which is the default.

# Returns
- A dictionary with six entries.
  - Key "1" => The first column (time).  
  - Key "2" => The second column (total mass - probability).  
  - Key "3" => The third column (SFR - original GADGET).  
  - Key "4" => The fourth column (SFR- probability).  
  - Key "5" => The fifth column (real total mass).  
  - Key "6" => The sixth column (real SFR).  
"""
function sfrTxtData(
    source_path::String,
    base_name::String;
    sim_cosmo::Int64 = 0,
    mass_unit::Unitful.FreeUnits = UnitfulAstro.Msun,
    time_unit::Unitful.FreeUnits = UnitfulAstro.Myr,
    sfr_unit::Unitful.FreeUnits = UnitfulAstro.Msun / UnitfulAstro.yr,
)::Dict{Union{Int64, String}, Any}

    # Get header of one snapshot for unit conversion.
    header = read_header(source_path * base_name * "_000")

    # Struct for unit conversion.
    if sim_cosmo == 1
        GU = GadgetPhysicalUnits(a_scale = header.time, hpar = header.h0)
        GU = @set GU.t_Myr = GU.t_Myr / sqrt(a)
        GU = @set GU.t_s = GU.t_s / sqrt(a)

        # Scale factor.
        a = header.time
        # Ω_K (curvature)
        omega_K = 1 - header.omega_0 - header.omega_l
        # Energy function.
        E = header.omega_0 / (a * a * a) + omega_K / (a * a) + header.omega_l
        # Hubble constant in 1 / Gyr.
        H = header.h0 * HUBBLE_CONST * sqrt(E)
        # Integrand of the time integral in Gyr.
        t_int = 1 / (H * a)
        # Conversion factor from internal units of time to physical units of time.
        t_conv = t_int * UnitfulAstro.Gyr
    else
        # For Newtonian simulation uses the default scale factor: a = 1.
        GU = GadgetPhysicalUnits(hpar = header.h0)

        # Conversion factor from internal units of time to physical units of time.
        t_conv = GU.t_Myr
    end

    # Get data from the file.
    sfr_txt = readdlm(source_path * "sfr.txt", Float64)

    # Column extraction and unit conversion.
    column_1 = @. ustrip(Float64, time_unit, sfr_txt[:, 1] * t_conv)
    column_2 = @. ustrip(Float64, mass_unit, sfr_txt[:, 2] * GU.m_cgs)
    column_3 = @. ustrip(
                    Float64, 
                    sfr_unit, 
                    sfr_txt[:, 3] * (UnitfulAstro.Msun / UnitfulAstro.yr)
                )
    column_4 = @. ustrip(
                    Float64, 
                    sfr_unit, 
                    sfr_txt[:, 4] * (UnitfulAstro.Msun / UnitfulAstro.yr)
                )
    column_5 = @. ustrip(Float64, mass_unit, sfr_txt[:, 5] * GU.m_cgs)
    column_6 = @. ustrip(Float64, sfr_unit, sfr_txt[:, 6] * (GU.m_cgs / GU.t_Myr))

    return Dict(
        1 => column_1,
        2 => column_2,
        3 => column_3,
        4 => column_4,
        5 => column_5,
        6 => column_6,
        "units" => Dict("mass" => mass_unit, "time" => time_unit, "sfr" => sfr_unit),
    )
end
