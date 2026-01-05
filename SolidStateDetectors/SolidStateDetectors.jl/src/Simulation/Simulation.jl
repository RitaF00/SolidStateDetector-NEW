abstract type AbstractSimulation{T<:SSDFloat} end

"""
    mutable struct Simulation{T <: SSDFloat, CS <: AbstractCoordinateSystem} <: AbstractSimulation{T}

Collection of all parts of a simulation of a [`SolidStateDetector`](@ref).

## Parametric types
* `T`: Precision type.
* `CS`: Coordinate system (`Cartesian` or `Cylindrical`).

## Fields 
* `config_dict::AbstractDict`: Dictionary (parsed configuration file) which initialized the simulation.
* `input_units::NamedTuple`: Units with which the `config_dict` should be parsed.
* `medium::NamedTuple`: Medium of the world.
* `detector::Union{SolidStateDetector{T}, Missing}`: The [`SolidStateDetector`](@ref) of the simulation.
* `world::World{T, 3, CS}`: The [`World`](@ref) of the simulation.
* `q_eff_imp::Union{EffectiveChargeDensity{T}, Missing}`: Effective charge resulting from the impurites in the [`Semiconductor`](@ref) of the `detector`.
* `imp_scale::Union{ImpurityScale{T}, Missing}`: Scale (alpha channel) of the impurity density (for depletion handling).  
* `q_eff_fix::Union{EffectiveChargeDensity{T}, Missing}`: Fixed charge resulting from fixed space charges in [`Passive`](@ref) of the `detector`.
* `ϵ_r::Union{DielectricDistribution{T}, Missing}`: The [`DielectricDistribution`](@ref) of the simulation.
* `point_types::Union{PointTypes{T}, Missing}`: The [`PointTypes`](@ref) of the simulation.
* `electric_potential::Union{ElectricPotential{T}, Missing}`: The [`ElectricPotential`](@ref) of the simulation.
* `weighting_potentials::Vector{Any}`: The [`WeightingPotential`](@ref) for each [`Contact`](@ref) of the `detector` in the simulation.
* `electric_field::Union{ElectricField{T}, Missing}`: The [`ElectricField`](@ref) of the simulation.
"""
mutable struct Simulation{T<:SSDFloat,CS<:AbstractCoordinateSystem} <: AbstractSimulation{T}
    config_dict::AbstractDict
    input_units::NamedTuple
    medium::NamedTuple # this should become a struct at some point
    detector::Union{SolidStateDetector{T},Missing}
    world::World{T,3,CS}
    q_eff_imp::Union{EffectiveChargeDensity{T},Missing} # Effective charge coming from the impurites of the semiconductors
    imp_scale::Union{ImpurityScale{T},Missing}
    q_eff_fix::Union{EffectiveChargeDensity{T},Missing} # Fixed charge coming from fixed space charges, e.g. charged up surface layers
    ϵ_r::Union{DielectricDistribution{T},Missing}
    point_types::Union{PointTypes{T},Missing}
    electric_potential::Union{ElectricPotential{T},Missing}
    weighting_potentials::Vector{Any}
    electric_field::Union{ElectricField{T},Missing}
end

function Simulation{T,CS}() where {T<:SSDFloat,CS<:AbstractCoordinateSystem}
    Simulation{T,CS}(
        Dict(),
        default_unit_tuple(),
        material_properties[materials["vacuum"]],
        missing,
        World(CS, (T(0), T(1), T(0), T(1), T(0), T(1))),
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        [missing],
        missing
    )
end

get_precision_type(::Simulation{T}) where {T} = T
get_coordinate_system(::Simulation{T,CS}) where {T,CS} = CS

function Base.NamedTuple(sim::Simulation{T}) where {T<:SSDFloat}
    wpots_strings = ["WeightingPotential_$(contact.id)" for contact in sim.detector.contacts]
    nt = (
        detector_json_string=_namedtuple(sim.config_dict),
        electric_potential=_namedtuple(sim.electric_potential),
        q_eff_imp=_namedtuple(sim.q_eff_imp),
        imp_scale=_namedtuple(sim.imp_scale),
        q_eff_fix=_namedtuple(sim.q_eff_fix),
        ϵ_r=_namedtuple(sim.ϵ_r),
        point_types=_namedtuple(sim.point_types),
        electric_field=_namedtuple(sim.electric_field),
        weighting_potentials=NamedTuple{Tuple(Symbol.(wpots_strings))}(_namedtuple.(sim.weighting_potentials))
    )
    return nt
end
Base.convert(T::Type{NamedTuple}, x::Simulation) = T(x)

function Simulation(nt::NamedTuple)
    missing_tuple = _namedtuple(missing)
    if nt.electric_potential !== missing_tuple
        epot = ElectricPotential(nt.electric_potential)
        T = eltype(epot.data)
        sim = Simulation{T}(_dict(nt.detector_json_string))
        sim.electric_potential = epot
        sim.q_eff_imp = EffectiveChargeDensity(nt.q_eff_imp)
        sim.q_eff_fix = EffectiveChargeDensity(nt.q_eff_fix)
        sim.ϵ_r = DielectricDistribution(nt.ϵ_r)
        sim.point_types = PointTypes(nt.point_types)
        sim.imp_scale = if !haskey(nt, :imp_scale)
            @warn """Stored simulation does not have a field for `imp_scale` (impurity scale) as this was 
            first introduced in SolidStateDetectors.jl v0.8 for improved depletion handling.
            It is advised to recalculate the simulation with the latest version.
            The field `imp_scale` is determined from `point_types`: 
                * undepleted point -> imp_scale = 0;  
                * depleted point -> imp_scale = 1;  
            """
            ImpurityScale(T.(.!is_undepleted_point_type.(sim.point_types.data)), sim.point_types.grid)
        else
            ImpurityScale(nt.imp_scale)
        end
        sim.electric_field = haskey(nt, :electric_field) && nt.electric_field !== missing_tuple ? ElectricField(nt.electric_field) : missing
    else
        T = Float32
        sim = Simulation{T}(_dict(nt.detector_json_string))
    end
    sim.weighting_potentials = if haskey(nt, :weighting_potentials)
        [
            let wp = Symbol("WeightingPotential_$(contact.id)")
                haskey(nt.weighting_potentials, wp) && getfield(nt.weighting_potentials, wp) !== missing_tuple ? WeightingPotential(getfield(nt.weighting_potentials, wp)) : missing
            end for contact in sim.detector.contacts
        ]
    else
        [missing for contact in sim.detector.contacts]
    end
    return sim
end
Base.convert(T::Type{Simulation}, x::NamedTuple) = T(x)

function Base.:(==)(sim1::P, sim2::P) where {P<:Union{Simulation,SolidStateDetector,AbstractObject,AbstractChargeDriftModel,AbstractTemperatureModel}}
    return typeof(sim1) == typeof(sim2) && all(broadcast(field -> isequal(getfield(sim1, field), getfield(sim2, field)), fieldnames(P)))
end


function println(io::IO, sim::Simulation{T}) where {T<:SSDFloat}
    println(typeof(sim), " - Coordinate system: ", get_coordinate_system(sim))
    println("  Environment Material: $(sim.medium.name)")
    println("  Detector: $(sim.detector.name)")
    println("  Electric potential: ", !ismissing(sim.electric_potential) ? size(sim.electric_potential) : missing)
    println("  Charge density: ", !ismissing(sim.q_eff_imp) ? size(sim.q_eff_imp) : missing)
    println("  Impurity scale: ", !ismissing(sim.imp_scale) ? size(sim.imp_scale) : missing)
    println("  Fix Charge density: ", !ismissing(sim.q_eff_fix) ? size(sim.q_eff_fix) : missing)
    println("  Dielectric distribution: ", !ismissing(sim.ϵ_r) ? size(sim.ϵ_r) : missing)
    println("  Point types: ", !ismissing(sim.point_types) ? size(sim.point_types) : missing)
    println("  Electric field: ", !ismissing(sim.electric_field) ? size(sim.electric_field) : missing)
    println("  Weighting potentials: ")
    for contact in sim.detector.contacts
        print("    Contact $(contact.id): ")
        println(!ismissing(sim.weighting_potentials[contact.id]) ? size(sim.weighting_potentials[contact.id]) : missing)
    end
end

function print(io::IO, sim::Simulation{T}) where {T<:SSDFloat}
    print(io, "Simulation{$T} - ", "$(sim.detector.name)")
end

function show(io::IO, sim::Simulation{T}) where {T<:SSDFloat}
    println(io, sim)
end

function show(io::IO, ::MIME"text/plain", sim::Simulation{T}) where {T<:SSDFloat}
    show(io, sim)
end


function Simulation{T}(dict::AbstractDict)::Simulation{T} where {T<:SSDFloat}
    CS::CoordinateSystemType = Cartesian
    if haskey(dict, "grid")
        if isa(dict["grid"], AbstractDict)
            CS = if dict["grid"]["coordinates"] == "cartesian"
                Cartesian
            elseif dict["grid"]["coordinates"] == "cylindrical"
                Cylindrical
            else
                @assert "`grid` in config file needs `coordinates` that are either `cartesian` or `cylindrical`"
            end
        elseif isa(dict["grid"], String)
            CS = if dict["grid"] == "cartesian"
                Cartesian
            elseif dict["grid"] == "cylindrical"
                Cylindrical
            else
                @assert "`grid` type in config file needs to be either `cartesian` or `cylindrical`"
            end
        end
    end
    sim::Simulation{T,CS} = Simulation{T,CS}()
    sim.config_dict = dict
    sim.input_units = construct_units(dict)
    sim.medium = material_properties[materials[haskey(dict, "medium") ? dict["medium"] : "vacuum"]]
    sim.detector = SolidStateDetector{T}(dict, sim.input_units)
    sim.world = if haskey(dict, "grid") && isa(dict["grid"], AbstractDict) && haskey(dict["grid"], "axes")
        World(T, dict["grid"], sim.input_units)
    else
        let det = sim.detector
            world_limits = get_world_limits_from_objects(CS, det)
            World(CS, world_limits)
        end
    end
    sim.weighting_potentials = Missing[missing for i in 1:length(sim.detector.contacts)]
    return sim
end

function Simulation{T}(config_file::AbstractString)::Simulation{T} where {T<:SSDFloat}
    dict = parse_config_file(config_file)
    return Simulation{T}(dict)
end
function Simulation(config_file::AbstractString)::Simulation{Float32}
    return Simulation{Float32}(config_file)
end

# Functions
"""
    Grid(sim::Simulation{T, Cartesian}; kwargs...)
    Grid(sim::Simulation{T, Cylindrical}; kwargs...)

Initializes a [`Grid`](@ref) based on the objects defined in a [`Simulation`](@ref).

The important points of all objects are sampled and added to the ticks of the grid.
The grid initialization can be tuned using a set of keyword arguments listed below.

## Arguments
* `sim::Simulation{T, S}`: [`Simulation`](@ref) for which the grid will be defined.

## Keywords
* `max_tick_distance = missing`: Maximum distance between neighbouring ticks of the grid.
    Additional grid ticks will be added if two neighbouring ticks are too far apart.
    `max_tick_distance` can either be a `Quantity`, e.g. `1u"mm"`, or a Tuple of `Quantity`, 
    e.g. `(1u"mm", 15u"°", 3u"mm")`,
    to set it for each axis of the `Grid` separately. Note that a `CartesianGrid3D` requires a 
    `Tuple{LengthQuantity, LengthQuantity, LengthQuantity}` while a `CylindricalGrid` requires a
    `Tuple{LengthQuantity, AngleQuantity, LengthQuantity}`.
    If `max_tick_distance` is `missing`, one fourth of the axis length is used.
* `max_distance_ratio::Real = 5`: If the ratio between a tick and its left and right neighbour
   is greater than `max_distance_ratio`, additional ticks are added between the ticks that are
   further apart. This prevents the ticks from being too unevenly spaced.
* `add_ticks_between_important_ticks::Bool = true`: If set to `true`, additional points
    will be added in between the important points obtained from sampling the objects of the
    simulation. If some objects are too close together, this will ensure a noticeable gap
    between them in the calculation of potentials and fields.
* `for_weighting_potential::Bool = false`: Grid will be optimized for the calculation of 
    an [`ElectricPotential`](@ref) if set to `true`, and of a [`WeightingPotential`](@ref)
    if set to `false`.
"""
function Grid(sim::Simulation{T,Cylindrical};
    for_weighting_potential::Bool=false,
    max_tick_distance::Union{Missing,LengthQuantity,Tuple{LengthQuantity,AngleQuantity,LengthQuantity}}=missing,
    max_distance_ratio::Real=5,
    add_ticks_between_important_ticks::Bool=true)::CylindricalGrid{T} where {T}
    det = sim.detector
    world = sim.world
    world_Δs = width.(world.intervals)
    world_Δr, world_Δφ, world_Δz = world_Δs

    samples::Vector{CylindricalPoint{T}} = sample(det, Cylindrical)
    important_r_ticks::Vector{T} = map(p -> p.r, samples)
    important_φ_ticks::Vector{T} = map(p -> p.φ, samples)
    important_z_ticks::Vector{T} = map(p -> p.z, samples)

    second_order_imp_ticks = if for_weighting_potential
        strong_electric_field_ticks = !ismissing(sim.electric_potential) ? get_ticks_at_positions_of_large_gradient(sim.electric_potential) : (T[], T[], T[])
        surface_of_depleted_volume_ticks = !ismissing(sim.imp_scale) ? get_ticks_at_positions_of_edge_of_depleted_volumes(sim.imp_scale) : (T[], T[], T[])
        vcat.(strong_electric_field_ticks, surface_of_depleted_volume_ticks)
    else
        (T[], T[], T[])
    end

    world_r_mid = mean(world.intervals[1])
    if for_weighting_potential && world_Δφ > 0
        world_φ_int = SSDInterval{T,:closed,:open,:periodic,:periodic}(0, 2π)
        world_Δφ = width(world_φ_int)
    else
        world_φ_int = world.intervals[2]
    end

    max_distance_z = T(world_Δz / 4)
    max_distance_φ = T(world_Δφ / 4)
    max_distance_r = T(world_Δr / 4)
    if !ismissing(max_tick_distance)
        if max_tick_distance isa LengthQuantity
            max_distance_z = max_distance_r = T(to_internal_units(max_tick_distance))
            max_distance_φ = max_distance_z / world_r_mid
        else #if max_tick_distance isa Tuple{LengthQuantity, AngleQuantity, LengthQuantity}
            max_distance_r = T(to_internal_units(max_tick_distance[1]))
            max_distance_φ = T(to_internal_units(max_tick_distance[2]))
            max_distance_z = T(to_internal_units(max_tick_distance[3]))
        end
    end

    println("Default maximum tick distance: z = $(max_distance_z), r = $(max_distance_r), ϕ = $(max_distance_φ) ")

    append!(important_r_ticks, endpoints(world.intervals[1])...)
    important_r_ticks = unique!(sort!(important_r_ticks))
    if add_ticks_between_important_ticks
        important_r_ticks = sort!(vcat(important_r_ticks, StatsBase.midpoints(important_r_ticks)))
    end
    iL = searchsortedfirst(important_r_ticks, world.intervals[1].left)
    iR = searchsortedfirst(important_r_ticks, world.intervals[1].right)
    important_r_ticks = unique(map(t -> isapprox(t, 0, atol=1e-12) ? zero(T) : t, important_r_ticks[iL:iR]))
    important_r_ticks = merge_close_ticks(important_r_ticks)
    imp2order_r_ticks = merge_close_ticks(second_order_imp_ticks[1], min_diff=world_Δs[1] / 20)
    important_r_ticks = merge_second_order_important_ticks(important_r_ticks, imp2order_r_ticks, min_diff=world_Δs[1] / 20)
    important_r_ticks = initialize_axis_ticks(important_r_ticks; max_ratio=T(max_distance_ratio))
    important_r_ticks = fill_up_ticks(important_r_ticks, max_distance_r)

    append!(important_z_ticks, endpoints(world.intervals[3])...)
    important_z_ticks = unique!(sort!(important_z_ticks))
    if add_ticks_between_important_ticks
        important_z_ticks = sort!(vcat(important_z_ticks, StatsBase.midpoints(important_z_ticks)))
    end
    iL = searchsortedfirst(important_z_ticks, world.intervals[3].left)
    iR = searchsortedfirst(important_z_ticks, world.intervals[3].right)
    important_z_ticks = unique(map(t -> isapprox(t, 0, atol=1e-12) ? zero(T) : t, important_z_ticks[iL:iR]))
    important_z_ticks = merge_close_ticks(important_z_ticks)
    imp2order_z_ticks = merge_close_ticks(second_order_imp_ticks[3], min_diff=world_Δs[3] / 20)
    important_z_ticks = merge_second_order_important_ticks(important_z_ticks, imp2order_z_ticks, min_diff=world_Δs[3] / 20)
    important_z_ticks = initialize_axis_ticks(important_z_ticks; max_ratio=T(max_distance_ratio))
    important_z_ticks = fill_up_ticks(important_z_ticks, max_distance_z)

    append!(important_φ_ticks, endpoints(world_φ_int)...)
    important_φ_ticks = unique!(sort!(important_φ_ticks))
    if add_ticks_between_important_ticks
        important_φ_ticks = sort!(vcat(important_φ_ticks, StatsBase.midpoints(important_φ_ticks)))
    end
    iL = searchsortedfirst(important_φ_ticks, world_φ_int.left)
    iR = searchsortedfirst(important_φ_ticks, world_φ_int.right)
    important_φ_ticks = unique(map(t -> isapprox(t, 0, atol=1e-3) ? zero(T) : t, important_φ_ticks[iL:iR]))
    important_φ_ticks = merge_close_ticks(important_φ_ticks, min_diff=T(1e-3))
    imp2order_φ_ticks = merge_close_ticks(second_order_imp_ticks[2], min_diff=world_Δs[2] / 20)
    important_φ_ticks = merge_second_order_important_ticks(important_φ_ticks, imp2order_φ_ticks, min_diff=world_Δs[2] / 20)
    important_φ_ticks = initialize_axis_ticks(important_φ_ticks; max_ratio=T(max_distance_ratio))
    important_φ_ticks = fill_up_ticks(important_φ_ticks, max_distance_φ)

    # r
    L, R, BL, BR = get_boundary_types(world.intervals[1])
    int_r = Interval{L,R,T}(endpoints(world.intervals[1])...)
    ax_r = even_tick_axis(DiscreteAxis{T,BL,BR}(int_r, important_r_ticks))

    # φ
    L, R, BL, BR = get_boundary_types(world_φ_int)
    int_φ = Interval{L,R,T}(endpoints(world_φ_int)...)
    ax_φ = if int_φ.left == int_φ.right
        DiscreteAxis{T,BL,BR}(int_φ, T[int_φ.left])
    else
        DiscreteAxis{T,BL,BR}(int_φ, important_φ_ticks)
    end
    if length(ax_φ) > 1
        φticks = if R == :open
            important_φ_ticks[1:end-1]
        else
            important_φ_ticks
        end
        ax_φ = typeof(ax_φ)(int_φ, φticks)
    end
    int_φ = ax_φ.interval
    if isodd(length(ax_φ)) && length(ax_φ) > 1 # must be even
        imax = findmax(diff(φticks))[2]
        push!(φticks, (φticks[imax] + φticks[imax+1]) / 2)
        sort!(φticks)
        ax_φ = typeof(ax_φ)(int_φ, φticks) # must be even
    end
    if length(ax_φ) > 1
        @assert iseven(length(ax_φ)) "CylindricalGrid must have even number of points in φ."
    end

    #z
    L, R, BL, BR = get_boundary_types(world.intervals[3])
    int_z = Interval{L,R,T}(endpoints(world.intervals[3])...)
    ax_z = even_tick_axis(DiscreteAxis{T,BL,BR}(int_z, important_z_ticks))

    total_points = length(ax_r) + length(ax_z)


    #println(" 📏​ The total number of points for each grid are: z = $(length(ax_z))")
    #println(" 📏​ The total number of points for each grid are: r = $(length(ax_r))")
    #println(" 📏​ The total number of points for each grid are: ϕ = $(length(ax_φ))")

    return CylindricalGrid{T}((ax_r, ax_φ, ax_z))
end


function Grid(sim::Simulation{T,Cartesian};
    max_tick_distance::Union{Missing,LengthQuantity,Tuple{LengthQuantity,LengthQuantity,LengthQuantity}}=missing,
    max_distance_ratio::Real=5,
    add_ticks_between_important_ticks::Bool=true,
    for_weighting_potential::Bool=false)::CartesianGrid3D{T} where {T}
    det = sim.detector
    world = sim.world
    world_Δs = width.(world.intervals)
    world_Δx, world_Δy, world_Δz = world_Δs

    samples::Vector{CartesianPoint{T}} = sample(det, Cartesian)
    important_x_ticks::Vector{T} = map(p -> p.x, samples)
    important_y_ticks::Vector{T} = map(p -> p.y, samples)
    important_z_ticks::Vector{T} = map(p -> p.z, samples)

    second_order_imp_ticks = if for_weighting_potential
        strong_electric_field_ticks = !ismissing(sim.electric_potential) ? get_ticks_at_positions_of_large_gradient(sim.electric_potential) : (T[], T[], T[])
        surface_of_depleted_volume_ticks = !ismissing(sim.imp_scale) ? get_ticks_at_positions_of_edge_of_depleted_volumes(sim.imp_scale) : (T[], T[], T[])
        vcat.(strong_electric_field_ticks, surface_of_depleted_volume_ticks)
    else
        (T[], T[], T[])
    end

    max_distance_x = T(world_Δx / 4)
    max_distance_y = T(world_Δy / 4)
    max_distance_z = T(world_Δz / 4)
    min_max_distance = min(max_distance_x, max_distance_y, max_distance_z)
    max_distance_x = max_distance_y = max_distance_z = min_max_distance
    if !ismissing(max_tick_distance)
        if max_tick_distance isa LengthQuantity
            max_distance_x = max_distance_y = max_distance_z =
                T(to_internal_units(max_tick_distance))
        else
            max_distance_x = T(to_internal_units(max_tick_distance[1]))
            max_distance_y = T(to_internal_units(max_tick_distance[2]))
            max_distance_z = T(to_internal_units(max_tick_distance[3]))
        end
    end

    append!(important_x_ticks, endpoints(world.intervals[1]))
    important_x_ticks = unique!(sort!(important_x_ticks))
    if add_ticks_between_important_ticks
        important_x_ticks = sort!(vcat(important_x_ticks, StatsBase.midpoints(important_x_ticks)))
    end
    iL = searchsortedfirst(important_x_ticks, world.intervals[1].left)
    iR = searchsortedfirst(important_x_ticks, world.intervals[1].right)
    important_x_ticks = unique(map(t -> isapprox(t, 0, atol=1e-12) ? zero(T) : t, important_x_ticks[iL:iR]))
    important_x_ticks = merge_close_ticks(important_x_ticks)
    imp2order_x_ticks = merge_close_ticks(second_order_imp_ticks[1], min_diff=world_Δs[1] / 20)
    important_x_ticks = merge_second_order_important_ticks(important_x_ticks, imp2order_x_ticks, min_diff=world_Δs[1] / 20)
    important_x_ticks = initialize_axis_ticks(important_x_ticks; max_ratio=T(max_distance_ratio))
    important_x_ticks = fill_up_ticks(important_x_ticks, max_distance_x)

    append!(important_y_ticks, endpoints(world.intervals[2]))
    important_y_ticks = unique!(sort!(important_y_ticks))
    if add_ticks_between_important_ticks
        important_y_ticks = sort!(vcat(important_y_ticks, StatsBase.midpoints(important_y_ticks)))
    end
    iL = searchsortedfirst(important_y_ticks, world.intervals[2].left)
    iR = searchsortedfirst(important_y_ticks, world.intervals[2].right)
    important_y_ticks = unique(map(t -> isapprox(t, 0, atol=1e-12) ? zero(T) : t, important_y_ticks[iL:iR]))
    important_y_ticks = merge_close_ticks(important_y_ticks)
    imp2order_y_ticks = merge_close_ticks(second_order_imp_ticks[2], min_diff=world_Δs[2] / 20)
    important_y_ticks = merge_second_order_important_ticks(important_y_ticks, imp2order_y_ticks, min_diff=world_Δs[2] / 20)
    important_y_ticks = initialize_axis_ticks(important_y_ticks; max_ratio=T(max_distance_ratio))
    important_y_ticks = fill_up_ticks(important_y_ticks, max_distance_y)

    append!(important_z_ticks, endpoints(world.intervals[3]))
    important_z_ticks = unique!(sort!(important_z_ticks))
    if add_ticks_between_important_ticks
        important_z_ticks = sort!(vcat(important_z_ticks, StatsBase.midpoints(important_z_ticks)))
    end
    iL = searchsortedfirst(important_z_ticks, world.intervals[3].left)
    iR = searchsortedfirst(important_z_ticks, world.intervals[3].right)
    important_z_ticks = unique(map(t -> isapprox(t, 0, atol=1e-12) ? zero(T) : t, important_z_ticks[iL:iR]))
    important_z_ticks = merge_close_ticks(important_z_ticks)
    imp2order_z_ticks = merge_close_ticks(second_order_imp_ticks[3], min_diff=world_Δs[3] / 20)
    important_z_ticks = merge_second_order_important_ticks(important_z_ticks, imp2order_z_ticks, min_diff=world_Δs[3] / 20)
    important_z_ticks = initialize_axis_ticks(important_z_ticks; max_ratio=T(max_distance_ratio))
    important_z_ticks = fill_up_ticks(important_z_ticks, max_distance_z)

    # x
    L, R, BL, BR = get_boundary_types(world.intervals[1])
    int_x = Interval{L,R,T}(endpoints(world.intervals[1])...)
    ax_x = even_tick_axis(DiscreteAxis{T,BL,BR}(int_x, important_x_ticks))

    # y
    L, R, BL, BR = get_boundary_types(world.intervals[2])
    int_y = Interval{L,R,T}(endpoints(world.intervals[2])...)
    ax_y = even_tick_axis(DiscreteAxis{T,BL,BR}(int_y, important_y_ticks))

    # z
    L, R, BL, BR = get_boundary_types(world.intervals[3])
    int_z = Interval{L,R,T}(endpoints(world.intervals[3])...)
    ax_z = even_tick_axis(DiscreteAxis{T,BL,BR}(int_z, important_z_ticks))

    return CartesianGrid3D{T}((ax_x, ax_y, ax_z))
end


function _guess_optimal_number_of_threads_for_SOR(gs::NTuple{3,Integer}, max_nthreads::Integer, S::Union{Type{Cylindrical},Type{Cartesian}})::Int
    max_nthreads = min(Base.Threads.nthreads(), max_nthreads)
    n = S == Cylindrical ? gs[2] * gs[3] : gs[1] * gs[2] # Number of grid points to be updated in each iteration of the outer loop
    return min(nextpow(2, max(cld(n + 1, 25), 4)), max_nthreads)
end


"""
    apply_initial_state!(sim::Simulation{T}, ::Type{ElectricPotential}, grid::Grid{T} = Grid(sim);
            not_only_paint_contacts::Bool = true, paint_contacts::Bool = true)::Nothing where {T <: SSDFloat}

Applies the initial state for the calculation of the [`ElectricPotential`](@ref).
It overwrites `sim.electric_potential`, `sim.q_eff_imp`, `sim.q_eff_fix`, `sim.ϵ` and `sim.point_types`
with the material properties and fixed potentials defined in `sim.detector`.

## Arguments 
* `sim::Simulation{T}`: [`Simulation`](@ref) for which the initial state should be applied.
* `grid::Grid{T}`: [`Grid`](@ref) to apply the initial state on. If no `grid` is given, 
    a default `Grid` is determined from `sim`.
    
## Keywords
* `not_only_paint_contacts::Bool = true`: Whether to only use the painting algorithm of the surfaces of [`Contact`](@ref)
    without checking if points are actually inside them.
    Setting it to `false` should improve the performance but the points inside of [`Contact`](@ref) are not fixed anymore.    
* `paint_contacts::Bool = true`: Enable or disable the painting of the surfaces of the [`Contact`](@ref) onto the `grid`.

## Examples
```julia
apply_initial_state!(sim, ElectricPotential, paint_contacts = false)
```
"""
function apply_initial_state!(sim::Simulation{T,CS}, ::Type{ElectricPotential}, grid::Grid{T}=Grid(sim);
    not_only_paint_contacts::Bool=true, paint_contacts::Bool=true)::Nothing where {T<:SSDFloat,CS}
    pcs = PotentialCalculationSetup(
        sim.detector, grid, sim.medium;
        use_nthreads=_guess_optimal_number_of_threads_for_SOR(size(grid), Base.Threads.nthreads(), CS),
        not_only_paint_contacts, paint_contacts
    )

    sim.q_eff_imp = EffectiveChargeDensity(EffectiveChargeDensityArray(pcs), grid)
    sim.imp_scale = ImpurityScale(ImpurityScaleArray(pcs), grid)
    sim.q_eff_fix = EffectiveChargeDensity(FixedEffectiveChargeDensityArray(pcs), grid)
    sim.ϵ_r = DielectricDistribution(DielectricDistributionArray(pcs), get_extended_midpoints_grid(grid))
    sim.point_types = PointTypes(PointTypeArray(pcs), grid)
    sim.electric_potential = ElectricPotential(ElectricPotentialArray(pcs), grid)
    nothing
end

"""
    apply_initial_state!(sim::Simulation{T}, ::Type{WeightingPotential}, contact_id::Int, grid::Grid{T} = Grid(sim))::Nothing

Applies the initial state for the calculation of the [`WeightingPotential`](@ref) for the [`Contact`}(@ref) with the id `contact_id`.
It overwrites `sim.weighting_potentials[contact_id]` with the fixed values on the [`Contact`}(@ref).

## Arguments 
* `sim::Simulation{T}`: [`Simulation`](@ref) for which the initial state should be applied.
* `contact_id::Int`: The `id` of the [`Contact`](@ref) for which the [`WeightingPotential`](@ref) is to be calculated.
* `grid::Grid{T}`: [`Grid`](@ref) to apply the initial state on. If no `grid` is given, 
    a default `Grid` is determined from `sim`.
    
## Keywords
* `not_only_paint_contacts::Bool = true`: Whether to only use the painting algorithm of the surfaces of [`Contact`](@ref)
    without checking if points are actually inside them.
    Setting it to `false` should improve the performance but the points inside of [`Contact`](@ref) are not fixed anymore.    
* `paint_contacts::Bool = true`: Enable or disable the painting of the surfaces of the [`Contact`](@ref) onto the `grid`.

## Examples
```julia
apply_initial_state!(sim, WeightingPotential, 1) # =>  applies initial state for weighting potential of contact with id 1
```
"""
function apply_initial_state!(sim::Simulation{T}, ::Type{WeightingPotential}, contact_id::Int, grid::Grid{T}=Grid(sim);
    not_only_paint_contacts::Bool=true, paint_contacts::Bool=true, depletion_handling::Bool=false)::Nothing where {T<:SSDFloat}
    pcs = PotentialCalculationSetup(
        sim.detector,
        grid,
        sim.medium,
        missing,
        depletion_handling ? sim.imp_scale.data : missing,
        weighting_potential_contact_id=contact_id;
        not_only_paint_contacts,
        paint_contacts,
        point_types=depletion_handling ? sim.point_types : missing
    )
    sim.weighting_potentials[contact_id] = WeightingPotential(ElectricPotentialArray(pcs), grid)
    nothing
end



"""
    update_till_convergence!( sim::Simulation{T} ::Type{ElectricPotential}, convergence_limit::Real; kwargs...)::T

Takes the current state of `sim.electric_potential` and updates it until it has converged.

There are several keyword arguments which can be used to tune the simulation.

## Arguments
* `sim::Simulation{T}`: [`Simulation`](@ref) for which `sim.electric_potential` will be updated.
* `convergence_limit::Real`: `convergence_limit` times the bias voltage sets the convergence limit of the relaxation.
    The convergence value is the absolute maximum difference of the potential between two iterations of all grid points.
    Default is `1e-7`.

## Keywords
* `n_iterations_between_checks::Int`: Number of iterations between checks. Default is set to `500`.
* `max_n_iterations::Int`: Set the maximum number of iterations which are performed after each grid refinement.
    Default is `-1`. If set to `-1` there will be no limit.
* `depletion_handling::Bool`: Enables the handling of undepleted regions. Default is `false`.
* `use_nthreads::Int`: Number of threads to use in the computation. Default is `Base.Threads.nthreads()`.
    The environment variable `JULIA_NUM_THREADS` must be set appropriately before the Julia session was
    started (e.g. `export JULIA_NUM_THREADS=8` in case of bash).
* `not_only_paint_contacts::Bool = true`: Whether to only use the painting algorithm of the surfaces of [`Contact`](@ref)
    without checking if points are actually inside them.
    Setting it to `false` should improve the performance but the points inside of [`Contact`](@ref) are not fixed anymore.    
* `paint_contacts::Bool = true`: Enable or disable the painting of the surfaces of the [`Contact`](@ref) onto the `grid`.
* `sor_consts::Union{<:Real, NTuple{2, <:Real}}`: Two element tuple in case of cylindrical coordinates.
    First element contains the SOR constant for `r` = 0.
    Second contains the constant at the outer most grid point in `r`. A linear scaling is applied in between.
    First element should be smaller than the second one and both should be `∈ [1.0, 2.0]`. Default is `[1.4, 1.85]`.
    In case of Cartesian coordinates, only one value is taken.
* `verbose::Bool=true`: Boolean whether info output is produced or not.
    
## Example 
```julia
SolidStateDetectors.update_till_convergence!(sim, ElectricPotential, 1e-6, depletion_handling = true)
```
"""
function update_till_convergence!(sim::Simulation{T,CS},
    ::Type{ElectricPotential},
    convergence_limit::Real=1e-7;
    n_iterations_between_checks::Int=500,
    max_n_iterations::Int=-1,
    depletion_handling::Bool=false,
    use_nthreads::Int=Base.Threads.nthreads(),
    not_only_paint_contacts::Bool=true,
    paint_contacts::Bool=true,
    device_array_type::Type{<:AbstractArray}=Array,
    sor_consts::Union{Missing,T,NTuple{2,T}}=missing,
    verbose::Bool=true
)::T where {T<:SSDFloat,CS<:AbstractCoordinateSystem}
    if ismissing(sor_consts)
        sor_consts = CS == Cylindrical ? (T(1.4), T(1.85)) : T(1.4)
    elseif length(sor_consts) == 1 && CS == Cylindrical
        sor_consts = (T(sor_consts), T(sor_consts))
    elseif length(sor_consts) > 1 && CS == Cartesian
        sor_consts = T(sor_consts[1])
    end
    only_2d = length(sim.electric_potential.grid.axes[2]) == 1

    pcs = Adapt.adapt(device_array_type, PotentialCalculationSetup(
        sim.detector, sim.electric_potential.grid, sim.medium, sim.electric_potential.data, sim.imp_scale.data, sor_consts=T.(sor_consts),
        use_nthreads=_guess_optimal_number_of_threads_for_SOR(size(sim.electric_potential.grid), Base.Threads.nthreads(), CS),
        not_only_paint_contacts=not_only_paint_contacts, paint_contacts=paint_contacts,
    ))

    via_KernelAbstractions = device_array_type <: GPUArrays.AnyGPUArray
    # This is just to be able to test the KernelAbstractions.jl backend on the CPU
    # as we cannot test it on GPU on GitHub. See also "SOR GPU Backend" test set.

    cf::T = _update_till_convergence!(pcs, T(convergence_limit), via_KernelAbstractions;
        only2d=Val{only_2d}(),
        depletion_handling=Val{depletion_handling}(),
        is_weighting_potential=Val{false}(),
        use_nthreads=use_nthreads,
        n_iterations_between_checks=n_iterations_between_checks,
        max_n_iterations=max_n_iterations,
        verbose=verbose)

    pcs = Adapt.adapt(Array, pcs)

    grid = Grid(pcs)
    sim.q_eff_imp = EffectiveChargeDensity(EffectiveChargeDensityArray(pcs), grid)
    sim.imp_scale = ImpurityScale(ImpurityScaleArray(pcs), grid)
    sim.q_eff_fix = EffectiveChargeDensity(FixedEffectiveChargeDensityArray(pcs), grid)
    sim.ϵ_r = DielectricDistribution(DielectricDistributionArray(pcs), get_extended_midpoints_grid(grid))
    sim.electric_potential = ElectricPotential(ElectricPotentialArray(pcs), grid)
    sim.point_types = PointTypes(PointTypeArray(pcs), grid)

    cf
end

"""
    update_till_convergence!( sim::Simulation{T} ::Type{WeightingPotential}, contact_id::Int, convergence_limit::Real; kwargs...)::T

Takes the current state of `sim.weighting_potentials[contact_id]` and updates it until it has converged.

There are several keyword arguments which can be used to tune the simulation.

## Arguments
* `sim::Simulation{T}`: [`Simulation`](@ref) for which `sim.weighting_potentials[contact_id]` will be updated.
* `contact_id::Int`: The `id` of the [`Contact`](@ref) for which the [`WeightingPotential`](@ref) is to be calculated.
* `convergence_limit::Real`: `convergence_limit` times the bias voltage sets the convergence limit of the relaxation.
    The convergence value is the absolute maximum difference of the potential between two iterations of all grid points.
    Default is `1e-7`.

## Keywords
* `n_iterations_between_checks::Int`: Number of iterations between checks. Default is set to `500`.
* `max_n_iterations::Int`: Set the maximum number of iterations which are performed after each grid refinement.
    Default is `-1`. If set to `-1` there will be no limit.
* `depletion_handling::Bool`: Enables the handling of undepleted regions. Default is `false`. This is an experimental feature:
    In undepleted regions (determined in `calculate_electric_potential!(sim; depletion_handling = true)`), the dielectric permittivity
    of the semiconductor is scaled up to mimic conductive behavior. The scale factor can be tuned via 
    the function [`scaling_factor_for_permittivity_in_undepleted_region`](@ref).
* `use_nthreads::Int`: Number of threads to use in the computation. Default is `Base.Threads.nthreads()`.
    The environment variable `JULIA_NUM_THREADS` must be set appropriately before the Julia session was
    started (e.g. `export JULIA_NUM_THREADS=8` in case of bash).
* `not_only_paint_contacts::Bool = true`: Whether to only use the painting algorithm of the surfaces of [`Contact`](@ref)
    without checking if points are actually inside them.
    Setting it to `false` should improve the performance but the points inside of [`Contact`](@ref) are not fixed anymore.    
* `paint_contacts::Bool = true`: Enable or disable the painting of the surfaces of the [`Contact`](@ref) onto the `grid`.
* `sor_consts::Union{<:Real, NTuple{2, <:Real}}`: Two element tuple in case of cylindrical coordinates.
    First element contains the SOR constant for `r` = 0.
    Second contains the constant at the outer most grid point in `r`. A linear scaling is applied in between.
    First element should be smaller than the second one and both should be `∈ [1.0, 2.0]`. Default is `[1.4, 1.85]`.
    In case of Cartesian coordinates, only one value is taken.
* `verbose::Bool=true`: Boolean whether info output is produced or not.
    
## Example 
```julia
SolidStateDetectors.update_till_convergence!(sim, WeightingPotential, 1, 1e-6, use_nthreads = 4)
```
"""
function update_till_convergence!(sim::Simulation{T,CS},
    ::Type{WeightingPotential},
    contact_id::Int,
    convergence_limit::Real=1e-7;
    n_iterations_between_checks::Int=500,
    max_n_iterations::Int=-1,
    depletion_handling::Bool=false,
    not_only_paint_contacts::Bool=true,
    paint_contacts::Bool=true,
    use_nthreads::Int=Base.Threads.nthreads(),
    device_array_type::Type{<:AbstractArray}=Array,
    sor_consts::Union{Missing,T,NTuple{2,T}}=missing,
    verbose::Bool=true
)::T where {T<:SSDFloat,CS<:AbstractCoordinateSystem}
    if ismissing(sor_consts)
        sor_consts = CS == Cylindrical ? (T(1.4), T(1.85)) : T(1.4)
    elseif length(sor_consts) == 1 && CS == Cylindrical
        sor_consts = (T(sor_consts), T(sor_consts))
    elseif length(sor_consts) > 1 && CS == Cartesian
        sor_consts = T(sor_consts[1])
    end

    only_2d::Bool = length(sim.weighting_potentials[contact_id].grid.axes[2]) == 1
    pcs = Adapt.adapt(device_array_type, PotentialCalculationSetup(
        sim.detector,
        sim.weighting_potentials[contact_id].grid,
        sim.medium,
        sim.weighting_potentials[contact_id].data,
        depletion_handling ? sim.imp_scale.data : ones(T, size(sim.weighting_potentials[contact_id].data)),
        sor_consts=T.(sor_consts),
        weighting_potential_contact_id=contact_id,
        use_nthreads=_guess_optimal_number_of_threads_for_SOR(size(sim.weighting_potentials[contact_id].grid), Base.Threads.nthreads(), CS),
        not_only_paint_contacts=not_only_paint_contacts,
        paint_contacts=paint_contacts,
        point_types=depletion_handling ? sim.point_types : missing)
    )

    via_KernelAbstractions = device_array_type <: GPUArrays.AnyGPUArray
    # This is just to be able to test the KernelAbstractions.jl backend on the CPU
    # as we cannot test it on GPU on GitHub. See also "SOR GPU Backend" test set.


    cf::T = _update_till_convergence!(pcs, T(convergence_limit), via_KernelAbstractions;
        only2d=Val{only_2d}(),
        depletion_handling=Val{false}(),
        is_weighting_potential=Val{true}(),
        use_nthreads=use_nthreads,
        n_iterations_between_checks=n_iterations_between_checks,
        max_n_iterations=max_n_iterations,
        verbose=verbose)

    pcs = Adapt.adapt(Array, pcs)
    sim.weighting_potentials[contact_id] = WeightingPotential(ElectricPotentialArray(pcs), sim.weighting_potentials[contact_id].grid)

    cf
end

"""
    refine!(sim::Simulation{T}, ::Type{ElectricPotential}, max_diffs::Tuple, minimum_distances::Tuple, kwargs...)

Takes the current state of `sim.electric_potential` and refines it with respect to the input arguments
`max_diffs` and `minimum_distances` by

1. extending the `grid` of `sim.electric_potential` to be a closed grid in all dimensions,
2. refining the axis of the grid based on `max_diffs` and `minimum_distances`:
   Insert new ticks between two existing ticks such that the potential difference between each tick becomes
   smaller than `max_diff[i]` (`i` -> dimension) but that the distances between the ticks stays larger than `minimum_distances[i]`, and
3. creating the new data array for the refined grid and fill it by interpolation of the the initial `grid`.


## Arguments
* `sim::Simulation{T}`: [`Simulation`](@ref) for which `sim.electric_potential` will be refined.
* `max_diffs::Tuple{<:Real,<:Real,<:Real}`: Maximum potential difference between two discrete ticks of `sim.electric_potential.grid` after refinement.
* `minimum_distances::Tuple{<:Real,<:Real,<:Real}`: Minimum distance (in SI Units) between two discrete ticks of `sim.electric_potential.grid` after refinement.
    
## Examples 
```julia 
SolidStateDetectors.refine!(sim, ElectricPotential, max_diffs = (100, 100, 100), minimum_distances = (0.01, 0.02, 0.01))
```
"""
function refine!(sim::Simulation{T}, ::Type{ElectricPotential},
    max_diffs::Tuple{<:Real,<:Real,<:Real}=(T(0), T(0), T(0)),
    minimum_distances::Tuple{<:Real,<:Real,<:Real}=(T(0), T(0), T(0));
    not_only_paint_contacts::Bool=true,
    paint_contacts::Bool=true,
    update_other_fields::Bool=false) where {T<:SSDFloat}
    sim.electric_potential = refine_scalar_potential(sim.electric_potential, T.(max_diffs), T.(minimum_distances))

    if update_other_fields
        pcs = PotentialCalculationSetup(sim.detector, sim.electric_potential.grid, sim.medium, sim.electric_potential.data,
            not_only_paint_contacts=not_only_paint_contacts, paint_contacts=paint_contacts)

        sim.imp_scale = ImpurityScale(ImpurityScaleArray(pcs), sim.electric_potential.grid)
        sim.q_eff_imp = EffectiveChargeDensity(EffectiveChargeDensityArray(pcs), sim.electric_potential.grid)
        sim.q_eff_fix = EffectiveChargeDensity(FixedEffectiveChargeDensityArray(pcs), sim.electric_potential.grid)
        sim.ϵ_r = DielectricDistribution(DielectricDistributionArray(pcs), get_extended_midpoints_grid(sim.electric_potential.grid))
        sim.point_types = PointTypes(PointTypeArray(pcs), sim.electric_potential.grid)
    end
    nothing
end


"""
    refine!(sim::Simulation{T}, ::Type{WeightingPotential}, max_diffs::Tuple{<:Real,<:Real,<:Real}, minimum_distances::Tuple{<:Real,<:Real,<:Real})

Takes the current state of `sim.weighting_potentials[contact_id]` and refines it with respect to the input arguments
`max_diffs` and `minimum_distances` by

1. extending the `grid` of `sim.weighting_potentials[contact_id]` to be a closed grid in all dimensions,
2. refining the axis of the grid based on `max_diffs` and `minimum_distances`:
   Insert new ticks between two existing ticks such that the potential difference between each tick becomes
   smaller than `max_diff[i]` (`i` -> dimension) but that the distances between the ticks stays larger than `minimum_distances[i]`, and
3. creating the new data array for the refined grid and fill it by interpolation of the the initial `grid`.


## Arguments
* `sim::Simulation{T}`: [`Simulation`](@ref) for which `sim.weighting_potentials[contact_id]` will be refined.
* `contact_id::Int`: The `id` of the [`Contact`](@ref) for which the [`WeightingPotential`](@ref) is refined.
* `max_diffs::Tuple{<:Real,<:Real,<:Real}`: Maximum potential difference between two discrete ticks of `sim.weighting_potentials[contact_id].grid` after refinement.
* `minimum_distances::Tuple{<:Real,<:Real,<:Real}`: Minimum distance (in SI Units) between two discrete ticks of `sim.weighting_potentials[contact_id].grid` after refinement.
   
## Examples 
```julia 
SolidStateDetectors.refine!(sim, WeightingPotential, 1, max_diffs = (0.01, 0.01, 0.01), minimum_distances = (0.01, 0.02, 0.01))
```
"""
function refine!(sim::Simulation{T}, ::Type{WeightingPotential}, contact_id::Int,
    max_diffs::Tuple{<:Real,<:Real,<:Real}=(T(0), T(0), T(0)),
    minimum_distances::Tuple{<:Real,<:Real,<:Real}=(T(0), T(0), T(0))) where {T<:SSDFloat}
    sim.weighting_potentials[contact_id] = refine_scalar_potential(sim.weighting_potentials[contact_id], max_diffs, minimum_distances)
    nothing
end

"""
    compute_min_tick_distance(grid)

Compute the minimum grid spacing (tick distance) for a grid.
Tick distances are computed as a small fraction of the physical axis lengths and then clamped to a safe numerical range.

# Arguments
- `grid` : A simulation grid (`CylindricalGrid` or `CartesianGrid3D`).

# Returns
A tuple of minimum tick distances in internal units, matching the grid type:
- Cylindrical: `(Δr, Δφ, Δz)`
- Cartesian: `(Δx, Δy, Δz)`
"""
function compute_min_tick_distance(grid::CylindricalGrid{T}) where {T}
    # qui vengono usati dei valori costanti
    min_tick = T(1e-12)
    max_tick = T(1e-5)
    fraction = T(1e-3)

    r_axis, phi_axis, z_axis = grid.axes

    r_len = width(r_axis.interval)
    z_len = width(z_axis.interval)

    safe_r_mid = max(abs(mean(r_axis.interval)), T(1e-30))

    Δr::T = clamp(r_len * fraction, min_tick, max_tick)
    Δφ::T = clamp(Δr / safe_r_mid, min_tick, max_tick)
    Δz::T = clamp(z_len * fraction, min_tick, max_tick)

    return (Δr, Δφ, Δz)
end

function compute_min_tick_distance(grid::CartesianGrid3D{T}) where {T}
    min_tick = T(1e-12)
    max_tick = T(1e-5)
    fraction = T(1e-3)

    x_axis, y_axis, z_axis = grid.axes

    x_len = width(x_axis.interval)
    y_len = width(y_axis.interval)
    z_len = width(z_axis.interval)

    Δx::T = clamp(x_len * fraction, min_tick, max_tick)
    Δy::T = clamp(y_len * fraction, min_tick, max_tick)
    Δz::T = clamp(z_len * fraction, min_tick, max_tick)

    return (Δx, Δy, Δz)
end


function _calculate_potential!(
    sim::Simulation{T,CS},
    potential_type::UnionAll,
    contact_id::Union{Missing,Int}=missing;
    convergence_limit::Real=1e-7,
    refinement_limits::Union{
        Missing,<:Real,Vector{<:Real},
        Tuple{<:Real,<:Real,<:Real},
        Vector{<:Tuple{<:Real,<:Real,<:Real}}
    }=[0.2, 0.1, 0.05],
    min_tick_distance::Union{
        Missing,LengthQuantity,
        Tuple{LengthQuantity,<:Union{LengthQuantity,AngleQuantity},LengthQuantity}
    }=missing,
    max_tick_distance::Union{
        Missing,LengthQuantity,
        Tuple{LengthQuantity,<:Union{LengthQuantity,AngleQuantity},LengthQuantity}
    }=missing,
    max_distance_ratio::Real=5,
    depletion_handling::Bool=false,
    use_nthreads::Union{Int,Vector{Int}}=Base.Threads.nthreads(),
    sor_consts::Union{Missing,<:Real,Tuple{<:Real,<:Real}}=missing,
    max_n_iterations::Int=50000,
    n_iterations_between_checks::Int=1000,
    not_only_paint_contacts::Bool=true,
    paint_contacts::Bool=true,
    verbose::Bool=true,
    device_array_type::Type{<:AbstractArray}=Array,
    initialize::Bool=true,
    grid::Union{Missing,Grid{T}}=initialize ? missing :
                                 (potential_type == ElectricPotential ?
                                  sim.electric_potential.grid :
                                  sim.weighting_potentials[contact_id].grid)
)::Nothing where {T<:SSDFloat,CS<:AbstractCoordinateSystem}

    begin # preparations
        onCPU = !(device_array_type <: GPUArrays.AnyGPUArray)
        convergence_limit = T(convergence_limit)
        isEP = potential_type == ElectricPotential
        isWP = !isEP

        if ismissing(grid)
            grid = Grid(sim;
                for_weighting_potential=isWP,
                max_tick_distance=max_tick_distance,
                max_distance_ratio=max_distance_ratio)
        end

        if ismissing(sor_consts)
            sor_consts = CS == Cylindrical ? (T(1.4), T(1.85)) : T(1.4)
        elseif length(sor_consts) == 1 && CS == Cylindrical
            sor_consts = (T(sor_consts), T(sor_consts))
        elseif length(sor_consts) > 1 && CS == Cartesian
            sor_consts = T(sor_consts[1])
        else
            sor_consts = T.(sor_consts)
        end

        new_min_tick_distance::NTuple{3,T} = begin
            if ismissing(min_tick_distance)
                compute_min_tick_distance(grid)
            elseif min_tick_distance isa LengthQuantity
                min_d = T(to_internal_units(min_tick_distance))
                if CS == Cylindrical
                    world_r_mid = mean(sim.world.intervals[1])
                    min_d, min_d / world_r_mid, min_d
                else
                    min_d, min_d, min_d
                end
            else
                T(to_internal_units(min_tick_distance[1])),
                T(to_internal_units(min_tick_distance[2])),
                T(to_internal_units(min_tick_distance[3]))
            end
        end

        refine = !ismissing(refinement_limits)
        refinement_limits isa Vector || (refinement_limits = [refinement_limits])
        n_refinement_steps = length(refinement_limits)
    end

    if initialize
        if isEP
            apply_initial_state!(sim, potential_type, grid;
                not_only_paint_contacts, paint_contacts)
            update_till_convergence!(sim, potential_type, convergence_limit;
                n_iterations_between_checks,
                max_n_iterations,
                depletion_handling,
                device_array_type,
                use_nthreads,
                sor_consts)
        else
            apply_initial_state!(sim, potential_type, contact_id, grid;
                not_only_paint_contacts, paint_contacts, depletion_handling)
            update_till_convergence!(sim, potential_type, contact_id, convergence_limit;
                n_iterations_between_checks,
                max_n_iterations,
                depletion_handling,
                device_array_type,
                use_nthreads,
                sor_consts)
        end
    end

    # ================= REFINEMENT LOOP =================
    if refine
        for iref in 1:n_refinement_steps
            is_last_ref = iref >= 3 || iref == n_refinement_steps
            ref_limits = T.(_extend_refinement_limits(refinement_limits[iref]))

            # ---------- GRID SIZE BEFORE ----------
            grid_before = isEP ?
                          sim.electric_potential.grid :
                          sim.weighting_potentials[contact_id].grid

            n_before = length(grid_before.axes[1]) +
                       length(grid_before.axes[3])

            # ---------- REFINE ----------
            if isEP
                bias_voltage = maximum(c.potential for c in sim.detector.contacts) -
                               minimum(c.potential for c in sim.detector.contacts)

                max_diffs = iszero(bias_voltage) ?
                            abs.(ref_limits .* (extrema(sim.electric_potential.data) |> x -> x[2] - x[1])) :
                            abs.(ref_limits .* bias_voltage)

                refine!(sim, ElectricPotential, max_diffs, new_min_tick_distance)
            else
                max_diffs = abs.(ref_limits)
                refine!(sim, WeightingPotential, contact_id, max_diffs, new_min_tick_distance)
            end

            # ---------- GRID SIZE AFTER ----------
            grid_after = isEP ?
                         sim.electric_potential.grid :
                         sim.weighting_potentials[contact_id].grid

            n_after = length(grid_after.axes[1]) +
                      length(grid_after.axes[3])

            println(
                "Refinement [$iref] → [$(iref+1)]: numero totale punti ",
                n_before == n_after ? "NON è cambiato" : "È CAMBIATO",
                " ($n_before → $n_after)"
            )


            # === BLOCCO DI GESTIONE POST-REFINE ===

            do_manual_wp_refine = (!isEP) && (iref <= 2)

            if n_before == n_after && do_manual_wp_refine

                #====
                # -------------------------------------------------
                # Raffinamento MANUALE asse z (solo WP, early steps)
                # -------------------------------------------------
                ax_z_old = grid_after.axes[3]
                z_ticks_old = ax_z_old.ticks

                println("📏 WP | punti z prima: ", length(z_ticks_old))

                # Trova i due intervalli più grandi
                diffs = diff(z_ticks_old)
                max_indices = partialsortperm(diffs, 1:2; rev=true)

                # Nuovi punti centrali
                new_z = [(z_ticks_old[i] + z_ticks_old[i+1]) / 2 for i in max_indices]
                println("➕ WP | nuovi punti z: ", new_z)

                # Nuovo asse z
                z_ticks_new = sort!(vcat(z_ticks_old, new_z))
                ax_z_new = typeof(ax_z_old)(ax_z_old.interval, z_ticks_new)

                println("📏 WP | punti z dopo: ", length(z_ticks_new))

                # Nuova griglia WP (solo asse z modificato)
                axes_new = (grid_after.axes[1], grid_after.axes[2], ax_z_new)
                grid_wp_new = Grid{eltype(grid_after),3,typeof(grid_after).parameters[3]}(axes_new)
                # ⚠️ Dati RI-ALLOCAZIONE sicura (NO riuso diretto)
                wp_old = sim.weighting_potentials[contact_id]
                wp_new_data = similar(wp_old.data, size(wp_old.data, 1),
                    size(wp_old.data, 2),
                    length(z_ticks_new))
                fill!(wp_new_data, zero(eltype(wp_new_data)))


                wp_new = WeightingPotential(wp_new_data, grid_wp_new)
                sim.weighting_potentials[contact_id] = wp_new

                println("✅ WP | asse z manuale applicato")

                # ===== PRINT DIAGNOSTICO (QUI) =====
                n_r = length(grid_wp_new.axes[1])
                n_phi = length(grid_wp_new.axes[2])
                n_z = length(grid_wp_new.axes[3])
                n_tot = n_r * n_phi * n_z

                println(
                    "🧲 WP | griglia MANUALE usata → ",
                    "r=$n_r, ϕ=$n_phi, z=$n_z (tot=$n_tot)"
                )

                # -------------------------------------------------
                # Raffinamento MANUALE asse r (solo WP, early steps)
                # -------------------------------------------------
                ax_r_old = grid_after.axes[1]
                r_ticks_old = ax_r_old.ticks

                println("📏 WP | punti r prima: ", length(r_ticks_old))

                # Trova i due intervalli più grandi
                diffs = diff(r_ticks_old)
                max_indices = partialsortperm(diffs, 1:2; rev=true)

                # Nuovi punti centrali
                new_r = [(r_ticks_old[i] + r_ticks_old[i+1]) / 2 for i in max_indices]
                println("➕ WP | nuovi punti r: ", new_r)

                # Nuovo asse r
                r_ticks_new = sort!(vcat(r_ticks_old, new_r))
                ax_r_new = typeof(ax_r_old)(ax_r_old.interval, r_ticks_new)

                println("📏 WP | punti r dopo: ", length(r_ticks_new))

                # Nuova griglia WP (solo asse r modificato)
                axes_new = (ax_r_new, grid_after.axes[2], grid_after.axes[3])
                grid_wp_new = Grid{eltype(grid_after),3,typeof(grid_after).parameters[3]}(axes_new)

                # ⚠️ Dati RI-ALLOCAZIONE sicura (NO riuso diretto)
                wp_old = sim.weighting_potentials[contact_id]
                wp_new_data = similar(
                    wp_old.data,
                    length(r_ticks_new),
                    size(wp_old.data, 2),
                    size(wp_old.data, 3)
                )
                fill!(wp_new_data, zero(eltype(wp_new_data)))

                wp_new = WeightingPotential(wp_new_data, grid_wp_new)
                sim.weighting_potentials[contact_id] = wp_new

                println("✅ WP | asse r manuale applicato")

                # ===== PRINT DIAGNOSTICO =====
                n_r = length(grid_wp_new.axes[1])
                n_phi = length(grid_wp_new.axes[2])
                n_z = length(grid_wp_new.axes[3])
                n_tot = n_r * n_phi * n_z

                println(
                    "🧲 WP | griglia MANUALE usata → ",
                    "r=$n_r, ϕ=$n_phi, z=$n_z (tot=$n_tot)"
                )

                ====================================#

                # -------------------------------------------------
                # Raffinamento MANUALE assi r + z (solo WP)
                # -------------------------------------------------

                # ===== ASSE r =====
                ax_r_old = grid_after.axes[1]
                r_ticks_old = ax_r_old.ticks
                println("📏 WP | punti r prima: ", length(r_ticks_old))

                diffs_r = diff(r_ticks_old)
                idx_r = partialsortperm(diffs_r, 1:2; rev=true)
                new_r = [(r_ticks_old[i] + r_ticks_old[i+1]) / 2 for i in idx_r]

                r_ticks_new = sort!(vcat(r_ticks_old, new_r))
                ax_r_new = typeof(ax_r_old)(ax_r_old.interval, r_ticks_new)

                println("➕ WP | nuovi punti r: ", new_r)
                println("📏 WP | punti r dopo: ", length(r_ticks_new))


                # ===== ASSE z =====
                ax_z_old = grid_after.axes[3]
                z_ticks_old = ax_z_old.ticks
                println("📏 WP | punti z prima: ", length(z_ticks_old))

                diffs_z = diff(z_ticks_old)
                idx_z = partialsortperm(diffs_z, 1:2; rev=true)
                new_z = [(z_ticks_old[i] + z_ticks_old[i+1]) / 2 for i in idx_z]

                z_ticks_new = sort!(vcat(z_ticks_old, new_z))
                ax_z_new = typeof(ax_z_old)(ax_z_old.interval, z_ticks_new)

                println("➕ WP | nuovi punti z: ", new_z)
                println("📏 WP | punti z dopo: ", length(z_ticks_new))


                # ===== NUOVA GRIGLIA WP =====
                axes_new = (ax_r_new, grid_after.axes[2], ax_z_new)
                grid_wp_new = Grid{eltype(grid_after),3,typeof(grid_after).parameters[3]}(axes_new)

                # ===== RI-ALLOCAZIONE DATI WP =====
                wp_old = sim.weighting_potentials[contact_id]

                wp_new_data = similar(
                    wp_old.data,
                    length(r_ticks_new),
                    size(wp_old.data, 2),
                    length(z_ticks_new)
                )
                fill!(wp_new_data, zero(eltype(wp_new_data)))

                wp_new = WeightingPotential(wp_new_data, grid_wp_new)
                sim.weighting_potentials[contact_id] = wp_new

                println("✅ WP | assi r + z manuali applicati")

                # ===== PRINT DIAGNOSTICO =====
                n_r = length(grid_wp_new.axes[1])
                n_phi = length(grid_wp_new.axes[2])
                n_z = length(grid_wp_new.axes[3])
                n_tot = n_r * n_phi * n_z

                println(
                    "🧲 WP | griglia MANUALE usata → ",
                    "r=$n_r, ϕ=$n_phi, z=$n_z (tot=$n_tot)"
                )

            else
                # -------------------------------------------------
                # refine! ha già cambiato la griglia → path standard
                # -------------------------------------------------
                println("🔁 Griglia raffinata automaticamente → nessun intervento manuale")
            end


            # da qui rinizia il vecchio codice 
            # ---------- SOLVE AGAIN ----------
            if isEP
                update_till_convergence!(sim, potential_type, convergence_limit;
                    n_iterations_between_checks,
                    max_n_iterations,
                    depletion_handling,
                    device_array_type,
                    use_nthreads,
                    sor_consts=is_last_ref ? T(1) : sor_consts)
            else
                update_till_convergence!(sim, potential_type, contact_id, convergence_limit;
                    n_iterations_between_checks,
                    max_n_iterations,
                    depletion_handling,
                    device_array_type,
                    use_nthreads,
                    sor_consts=is_last_ref ? T(1) : sor_consts)
            end
        end
    end

    nothing
end



"""
    calculate_weighting_potential!(sim::Simulation{T}, contact_id::Int; kwargs...)::Nothing

Calculates the [`WeightingPotential`](@ref) for a [`Contact`](@ref) with `contact_id` 
given [`Simulation`](@ref) `sim` on an adaptive grid through successive over relaxation 
and stores it in `sim.weighting_potentials[contact_id]`.

There are several keyword arguments which can be used to tune the calculation.

## Keywords
* `convergence_limit::Real`: `convergence_limit` sets the convergence limit of the relaxation.
    The convergence value is the absolute maximum difference of the potential between two iterations of all grid points.
    Default of `convergence_limit` is `1e-7`.
* `refinement_limits`: Defines the maximum relative allowed differences 
    of the potential value of neighbored grid points 
    in each dimension for each refinement.
    - `rl::Real` -> One refinement with `rl` equal in all 3 dimensions.
    - `rl::Tuple{<:Real,<:Real,<:Real}` -> One refinement with `rl` set individual for each dimension.
    - `rl::Vector{<:Real}` -> `length(l)` refinements with `rl[i]` being the limit for the i-th refinement. 
    - `rl::Vector{<:Real,<:Real,<:Real}}` -> `length(rl)` refinements with `rl[i]` being the limits for the `i`-th refinement.
* `min_tick_distance::Tuple{<:Quantity, <:Quantity, <:Quantity}`: Tuple of the minimum allowed distance between 
    two grid ticks for each dimension. It prevents the refinement to make the grid too fine.
    Default is `1e-5` for linear axes and `1e-5 / (0.25 * r_max)` for the polar axis in case of a cylindrical `grid`.
* `max_tick_distance::Tuple{<:Quantity, <:Quantity, <:Quantity}`: Tuple of the maximum allowed distance between 
    two grid ticks for each dimension used in the initialization of the grid.
    Default is 1/4 of size of the world of the respective dimension.
* `max_distance_ratio::Real`: Maximum allowed ratio between the two distances in any dimension to the two neighbouring grid points. 
        If the ratio is too large, additional ticks are generated such that the new ratios are smaller than `max_distance_ratio`.
        Default is `5`.
* `grid::Grid`: Initial grid used to start the simulation. Default is `Grid(sim)`.
* `depletion_handling::Bool`: Enables the handling of undepleted regions. Default is `false`. This is an experimental feature:
    In undepleted regions (determined in `calculate_electric_potential!(sim; depletion_handling = true)`), the dielectric permittivity
    of the semiconductor is scaled up to mimic conductive behavior. The scale factor can be tuned via 
    the function [`scaling_factor_for_permittivity_in_undepleted_region`](@ref).
* `use_nthreads::Union{Int, Vector{Int}}`: If `<:Int`, `use_nthreads` defines the maximum number of threads to be used in the computation. 
    Fewer threads might be used depending on the current grid size due to threading overhead. Default is `Base.Threads.nthreads()`.
    If `<:Vector{Int}`, `use_nthreads[i]` defines the number of threads used for each grid (refinement) stage of the field simulation.
    The environment variable `JULIA_NUM_THREADS` must be set appropriately before the Julia session was
    started (e.g. `export JULIA_NUM_THREADS=8` in case of bash).
* `sor_consts::Union{<:Real, NTuple{2, <:Real}}`: Two element tuple in case of cylindrical coordinates.
    First element contains the SOR constant for `r` = 0.
    Second contains the constant at the outer most grid point in `r`. A linear scaling is applied in between.
    First element should be smaller than the second one and both should be `∈ [1.0, 2.0]`. Default is `[1.4, 1.85]`.
    In case of Cartesian coordinates, only one value is taken.
* `max_n_iterations::Int`: Set the maximum number of iterations which are performed after each grid refinement.
    Default is `50000`. If set to `-1` there will be no limit.
* `not_only_paint_contacts::Bool = true`: Whether to only use the painting algorithm of the surfaces of [`Contact`](@ref)
    without checking if points are actually inside them.
    Setting it to `false` should improve the performance but the points inside of [`Contact`](@ref) are not fixed anymore.    
* `paint_contacts::Bool = true`: Enable or disable the painting of the surfaces of the [`Contact`](@ref) onto the `grid`.
* `verbose::Bool=true`: Boolean whether info output is produced or not.

## Example 
```julia 
calculate_weighting_potential!(sim, 1, refinement_limits = [0.3, 0.1, 0.05], max_distance_ratio = 4, max_n_iterations = 20000)
```
"""
function calculate_weighting_potential!(sim::Simulation{T}, contact_id::Int, args...; n_points_in_φ::Union{Missing,Int}=missing, kwargs...)::Nothing where {T<:SSDFloat}
    _calculate_potential!(sim, WeightingPotential, contact_id, args...; kwargs...)
    nothing
end


"""
    calculate_electric_potential!(sim::Simulation{T}; kwargs...)::Nothing


Calculates the [`ElectricPotential`](@ref) for a given [`Simulation`](@ref) `sim` on an adaptive grid
through successive over relaxation and stores it in `sim.electric_potential`.

There are several keyword arguments which can be used to tune the calculation.

## Arguments
* `sim::Simulation{T}`: [`Simulation`](@ref) for which the [`ElectricPotential`](@ref) is calculated.

## Keywords
* `convergence_limit::Real`: `convergence_limit` times the bias voltage sets the convergence limit of the relaxation.
    The convergence value is the absolute maximum difference of the potential between two iterations of all grid points.
    Default of `convergence_limit` is `1e-7` (times bias voltage).
* `refinement_limits`: Defines the maximum relative (to applied bias voltage) allowed differences 
    of the potential value of neighbored grid points 
    in each dimension for each refinement.
    - `rl::Real` -> One refinement with `rl` equal in all 3 dimensions.
    - `rl::Tuple{<:Real,<:Real,<:Real}` -> One refinement with `rl` set individual for each dimension.
    - `rl::Vector{<:Real}` -> `length(l)` refinements with `rl[i]` being the limit for the i-th refinement. 
    - `rl::Vector{<:Real,<:Real,<:Real}}` -> `length(rl)` refinements with `rl[i]` being the limits for the `i`-th refinement.
* `min_tick_distance::Tuple{<:Quantity, <:Quantity, <:Quantity}`: Tuple of the minimum allowed distance between 
    two grid ticks for each dimension. It prevents the refinement to make the grid too fine.
    Default is `1e-5` for linear axes and `1e-5 / (0.25 * r_max)` for the polar axis in case of a cylindrical `grid`.
* `max_tick_distance::Tuple{<:Quantity, <:Quantity, <:Quantity}`: Tuple of the maximum allowed distance between 
    two grid ticks for each dimension used in the initialization of the grid.
    Default is 1/4 of size of the world of the respective dimension.
* `max_distance_ratio::Real`: Maximum allowed ratio between the two distances in any dimension to the two neighbouring grid points. 
        If the ratio is too large, additional ticks are generated such that the new ratios are smaller than `max_distance_ratio`.
        Default is `5`.
* `grid::Grid`: Initial grid used to start the simulation. Default is `Grid(sim)`.
* `depletion_handling::Bool`: Enables the handling of undepleted regions. Default is `false`.
* `use_nthreads::Union{Int, Vector{Int}}`: If `<:Int`, `use_nthreads` defines the maximum number of threads to be used in the computation. 
    Fewer threads might be used depending on the current grid size due to threading overhead. Default is `Base.Threads.nthreads()`.
    If `<:Vector{Int}`, `use_nthreads[i]` defines the number of threads used for each grid (refinement) stage of the field simulation.
    The environment variable `JULIA_NUM_THREADS` must be set appropriately before the Julia session was
    started (e.g. `export JULIA_NUM_THREADS=8` in case of bash).
* `sor_consts::Union{<:Real, NTuple{2, <:Real}}`: Two element tuple in case of cylindrical coordinates.
    First element contains the SOR constant for `r` = 0.
    Second contains the constant at the outer most grid point in `r`. A linear scaling is applied in between.
    First element should be smaller than the second one and both should be `∈ [1.0, 2.0]`. Default is `[1.4, 1.85]`.
    In case of Cartesian coordinates, only one value is taken.
* `max_n_iterations::Int`: Set the maximum number of iterations which are performed after each grid refinement.
    Default is `50000`. If set to `-1` there will be no limit.
* `not_only_paint_contacts::Bool = true`: Whether to only use the painting algorithm of the surfaces of [`Contact`](@ref)
    without checking if points are actually inside them.
    Setting it to `false` should improve the performance but the points inside of [`Contact`](@ref) are not fixed anymore.    
* `paint_contacts::Bool = true`: Enable or disable the painting of the surfaces of the [`Contact`](@ref) onto the `grid`.
* `verbose::Bool=true`: Boolean whether info output is produced or not.

## Example 
```julia 
calculate_electric_potential!(sim, refinement_limits = [0.3, 0.1, 0.05], max_distance_ratio = 4, max_n_iterations = 20000)
```
"""
function calculate_electric_potential!(sim::Simulation{T}, args...; kwargs...)::Nothing where {T<:SSDFloat}
    _calculate_potential!(sim, ElectricPotential, args...; kwargs...)
    nothing
end

"""
    calculate_electric_field!(sim::Simulation{T}; n_points_in_φ::Union{Missing, Int} = missing)::Nothing

Calculates the [`ElectricField`](@ref) from the [`ElectricPotential`](@ref) stored in `sim.electric_potential` and stores it in
`sim.electric_field`. 

## Arguments
* `sim::Simulation{T}`: [`Simulation`](@ref) for which `sim.electric_potential` has already been calculated.

## Keywords
* `n_points_in_φ::Union{Missing, Int}`: For a 2D [`ElectricPotential`](@ref) (cylindrical coordinates and symmetric in `φ`), `sim.electric_potential`
    is extended to `n_points_in_φ` "layers" in `φ` in order to calculate a 3D [`ElectricField`]. If `n_points_in_φ` is `missing`, the 
    default value is `36`.

## Examples 
    calculate_electric_field!(sim, n_points_in_φ = 32)

!!! note 
    This method only works if `sim.electric_potential` has already been calculated and is not `missing`.
"""
function calculate_electric_field!(sim::Simulation{T,CS}; n_points_in_φ::Union{Missing,Int}=missing, use_nthreads::Int=Base.Threads.nthreads())::Nothing where {T<:SSDFloat,CS}
    @assert !ismissing(sim.electric_potential) "Electric potential has not been calculated yet. Please run `calculate_electric_potential!(sim)` first."
    periodicity::T = width(sim.world.intervals[2])
    e_pot, point_types = if CS == Cylindrical && periodicity == T(0) # 2D, only one point in φ
        if ismissing(n_points_in_φ)
            @info "\tIn electric field calculation: Keyword `n_points_in_φ` not set.\n\t\tDefault is `n_points_in_φ = 36`. 2D field will be extended to 36 points in φ."
            n_points_in_φ = 36
        else
            if !(n_points_in_φ > 1 && iseven(n_points_in_φ))
                @info "\tIn electric field calculation: Keyword `n_points_in_φ` is $(n_points_in_φ) but must be even and larger than 1.\n\t\t`n_points_in_φ` is now set to 36. 2D field will be extended to 36 points in φ."
                n_points_in_φ = 36
            end
        end
        get_2π_potential(sim.electric_potential, n_points_in_φ=n_points_in_φ),
        get_2π_potential(sim.point_types, n_points_in_φ=n_points_in_φ)
    elseif CS == Cylindrical
        get_2π_potential(sim.electric_potential),
        get_2π_potential(sim.point_types)
    else
        sim.electric_potential,
        sim.point_types
    end
    sim.electric_field = get_electric_field_from_potential(e_pot, point_types; use_nthreads)
    nothing
end

function drift_charges(sim::Simulation{T}, starting_positions::VectorOfArrays{CartesianPoint{T}}, energies::VectorOfArrays{T};
    Δt::RealQuantity=5u"ns", max_nsteps::Int=1000, diffusion::Bool=false, self_repulsion::Bool=false,
    end_drift_when_no_field::Bool=true, geometry_check::Bool=false, verbose::Bool=true)::Vector{EHDriftPath{T}} where {T<:SSDFloat}
    return _drift_charges(sim.detector, sim.point_types.grid, sim.point_types, starting_positions, energies,
        interpolated_vectorfield(sim.electric_field), Δt;
        max_nsteps, diffusion, self_repulsion, end_drift_when_no_field, geometry_check, verbose)
end

function get_signal(sim::Simulation{T,CS}, drift_paths::Vector{EHDriftPath{T}}, energy_depositions::Vector{T}, contact_id::Int; Δt::TT=T(5) * u"ns", signal_unit::Unitful.Units=u"e_au") where {T<:SSDFloat,CS,TT}
    dt::T = to_internal_units(Δt)
    wpot::Interpolations.Extrapolation{T,3} = interpolated_scalarfield(sim.weighting_potentials[contact_id])
    timestamps = _common_timestamps(drift_paths, dt)
    signal::Vector{T} = zeros(T, length(timestamps))
    add_signal!(signal, timestamps, drift_paths, energy_depositions, wpot, sim.point_types, sim.detector.semiconductor.charge_trapping_model)
    calibration_factor::Quantity{T,dimension(signal_unit)} = _convert_internal_energy_to_external_unit(signal_unit, sim.detector.semiconductor.material)
    return RDWaveform(range(zero(T) * unit(Δt), step=T(ustrip(Δt)) * unit(Δt), length=length(signal)), signal * calibration_factor)
end

"""
    simulate!( sim::Simulation{T}; kwargs...) where {T, S}


Performs a full chain simulation for a given [`Simulation`](@ref) by
    
1. calculating the [`ElectricPotential`](@ref),
2. calculating the [`ElectricField`](@ref),
3. calculating the [`WeightingPotential`](@ref) for each [`Contact`](@ref).

The output is stored in `sim.electric_potential`, `sim.electric_field` and `sim.weighting_potentials`, respectively.

There are several keyword arguments which can be used to tune the simulation.

## Arguments
* `sim::Simulation{T}`: [`Simulation`](@ref) for which the full chain simulation should be performed.


## Keywords
* `convergence_limit::Real`: `convergence_limit` times the bias voltage sets the convergence limit of the relaxation.
    The convergence value is the absolute maximum difference of the potential between two iterations of all grid points.
    Default of `convergence_limit` is `1e-7` (times bias voltage).
* `refinement_limits`: Defines the maximum relative (to applied bias voltage) allowed differences 
    of the potential value of neighboured grid points 
    in each dimension for each refinement.
    - `rl::Real` -> One refinement with `rl` equal in all 3 dimensions.
    - `rl::Tuple{<:Real,<:Real,<:Real}` -> One refinement with `rl` set individual for each dimension.
    - `rl::Vector{<:Real}` -> `length(l)` refinements with `rl[i]` being the limit for the i-th refinement. 
    - `rl::Vector{<:Real,<:Real,<:Real}}` -> `length(rl)` refinements with `rl[i]` being the limits for the `i`-th refinement.
* `min_tick_distance::Tuple{<:Quantity, <:Quantity, <:Quantity}`: Tuple of the minimum allowed distance between 
    two grid ticks for each dimension. It prevents the refinement to make the grid too fine.
    Default is `1e-5` for linear axes and `1e-5 / (0.25 * r_max)` for the polar axis in case of a cylindrical `grid`.
* `max_tick_distance::Tuple{<:Quantity, <:Quantity, <:Quantity}`: Tuple of the maximum allowed distance between 
    two grid ticks for each dimension used in the initialization of the grid.
    Default is 1/4 of size of the world of the respective dimension.
* `max_distance_ratio::Real`: Maximum allowed ratio between the two distances in any dimension to the two neighbouring grid points. 
        If the ratio is too large, additional ticks are generated such that the new ratios are smaller than `max_distance_ratio`.
        Default is `5`.
* `depletion_handling::Bool`: Enables the handling of undepleted regions. Default is `false`.
* `use_nthreads::Int`: Number of threads to use in the computation. Default is `Base.Threads.nthreads()`.
    The environment variable `JULIA_NUM_THREADS` must be set appropriately before the Julia session was
    started (e.g. `export JULIA_NUM_THREADS=8` in case of bash).
* `sor_consts::Union{<:Real, NTuple{2, <:Real}}`: Two element tuple in case of cylindrical coordinates.
    First element contains the SOR constant for `r` = 0.
    Second contains the constant at the outer most grid point in `r`. A linear scaling is applied in between.
    First element should be smaller than the second one and both should be `∈ [1.0, 2.0]`. Default is `[1.4, 1.85]`.
    In case of Cartesian coordinates, only one value is taken.
* `max_n_iterations::Int`: Set the maximum number of iterations which are performed after each grid refinement.
    If set to `-1` there will be no limit. Default is no limit.
* `not_only_paint_contacts::Bool = true`: Whether to only use the painting algorithm of the surfaces of [`Contact`](@ref)
    without checking if points are actually inside them.
    Setting it to `false` should improve the performance but the points inside of [`Contact`](@ref) are not fixed anymore.    
* `paint_contacts::Bool = true`: Enable or disable the painting of the surfaces of the [`Contact`](@ref) onto the `grid`.
* `verbose::Bool=true`: Boolean whether info output is produced or not.

See also [`calculate_electric_potential!`](@ref), [`calculate_electric_field!`](@ref) and [`calculate_weighting_potential!`](@ref).

## Example 
```julia 
simulate!(sim, refinement_limits = [0.3, 0.1, 0.05], max_distance_ratio = 4, max_n_iterations = 20000)
```
"""
function simulate!(sim::Simulation{T,S};
    convergence_limit::Real=1e-7,
    refinement_limits=[0.2, 0.1, 0.05],
    min_tick_distance::Union{Missing,LengthQuantity,Tuple{LengthQuantity,AngleQuantity,LengthQuantity}}=missing,
    max_tick_distance::Union{Missing,LengthQuantity,Tuple{LengthQuantity,AngleQuantity,LengthQuantity}}=missing,
    max_distance_ratio::Real=5,
    depletion_handling::Bool=false,
    use_nthreads::Union{Int,Vector{Int}}=Base.Threads.nthreads(),
    sor_consts::Union{Missing,<:Real,Tuple{<:Real,<:Real}}=missing,
    max_n_iterations::Int=-1,
    device_array_type::Type{<:AbstractArray}=Array,
    not_only_paint_contacts::Bool=true,
    paint_contacts::Bool=true,
    verbose::Bool=false) where {T<:SSDFloat,S}
    calculate_electric_potential!(sim,
        convergence_limit=convergence_limit,
        refinement_limits=refinement_limits,
        min_tick_distance=min_tick_distance,
        max_tick_distance=max_tick_distance,
        max_distance_ratio=max_distance_ratio,
        depletion_handling=depletion_handling,
        use_nthreads=use_nthreads,
        sor_consts=sor_consts,
        max_n_iterations=max_n_iterations,
        device_array_type=device_array_type,
        not_only_paint_contacts=not_only_paint_contacts,
        paint_contacts=paint_contacts,
        verbose=verbose
    )
    for contact in sim.detector.contacts
        calculate_weighting_potential!(sim, contact.id,
            convergence_limit=convergence_limit,
            refinement_limits=refinement_limits,
            min_tick_distance=min_tick_distance,
            max_tick_distance=max_tick_distance,
            max_distance_ratio=max_distance_ratio,
            depletion_handling=depletion_handling,
            use_nthreads=use_nthreads,
            sor_consts=sor_consts,
            max_n_iterations=max_n_iterations,
            device_array_type=device_array_type,
            not_only_paint_contacts=not_only_paint_contacts,
            paint_contacts=paint_contacts,
            verbose=verbose
        )
    end
    calculate_electric_field!(sim)
    @info "Detector simulation done"
end

include("ElectricFieldEnergy.jl")
include("Capacitance.jl")
include("DepletionVoltage.jl")
