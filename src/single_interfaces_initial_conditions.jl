abstract type AbstractSingleInterfaceInitialConditions <: AbstractInitialConditions end

Base.iterate(sics::AbstractSingleInterfaceInitialConditions, state = 1) =
    state > length(fieldnames(sics)) ? nothing : (getfield(sics, state), state + 1)
"""
    struct STSingleInterfaceInitialConditions
Initial conditions for a single interface (i.e. two layers with uniform `S` and `T` seperated
by a step change). The property `maintain_interface` is a `Boolean` which if set to `true` will
set [reentrant_boundary_conditions](@ref) so that the interface will be maintained (by not
letting the system run down).
"""
struct STSingleInterfaceInitialConditions{T, A, IS} <: AbstractSingleInterfaceInitialConditions
    "The depth of the interface"
     depth_of_interface :: T
    "Salinity values in each layer"
        salinity_values :: A
     "Temperature values in each layer"
     temperature_values :: A
     "Initial smoothing function over interface"
    interface_smoothing :: IS
     "Initial R_ρ at the interface"
                    R_ρ :: T
end
function STSingleInterfaceInitialConditions(model, depth_of_interface, salinity, temperature;
                                            interface_smoothing = nothing)

    eos = model.buoyancy.model.equation_of_state

    R_ρ = compute_R_ρ(salinity, temperature, depth_of_interface, eos)

    return STSingleInterfaceInitialConditions(depth_of_interface, salinity, temperature,
                                              interface_smoothing, R_ρ)

end
function STSingleInterfaceInitialConditions(eos::BoussinesqEquationOfState, depth_of_interface,
                                            salinity, temperature; interface_smoothing = nothing)

    R_ρ = compute_R_ρ(salinity, temperature, depth_of_interface, eos)

    return STSingleInterfaceInitialConditions(depth_of_interface, salinity, temperature,
                                              interface_smoothing, R_ρ)

end
const SingleInterfaceICs = STSingleInterfaceInitialConditions # alias
Base.summary(ics::SingleInterfaceICs) = "Single S-T interface at z = $(ics.depth_of_interface) on bounded z domain"

"""
    struct PeriodicSTSingleInterfaceInitialConditions
Sets a `BackgroundField` according to `background_State` and uses a triply periodic domain
to evolve salinity and temperature anomalies about the background state.
"""
struct PeriodicSTSingleInterfaceInitialConditions{T, A, IS, BF} <: AbstractSingleInterfaceInitialConditions
    "The depth of the interface"
     depth_of_interface :: T
    "Salinity values in each layer"
        salinity_values :: A
     "Temperature values in each layer"
     temperature_values :: A
    "Initial smoothing function over interface"
    interface_smoothing :: IS
     "Initial R_ρ at the interface"
                    R_ρ :: T
    "`BackgroundField` about which anomaly is advected. Should be an `AbstractBackgroundFunction`."
       background_state :: BF
end
function PeriodicSTSingleInterfaceInitialConditions(eos::BoussinesqEquationOfState, depth_of_interface,
                                                    salinity, temperature;
                                                    interface_smoothing = NoSmoothing,
                                                    background_state = NoBackground)

    R_ρ = compute_R_ρ(salinity, temperature, depth_of_interface, eos)

    return PeriodicSTSingleInterfaceInitialConditions(depth_of_interface, salinity, temperature,
                                                      interface_smoothing, R_ρ, background_state)

end
const PeriodoicSingleInterfaceICs = PeriodicSTSingleInterfaceInitialConditions # alias
Base.summary(ics::PeriodoicSingleInterfaceICs) = "Single S-T interface at z = $(ics.depth_of_interface) on triply periodic domain"

function compute_R_ρ(salinity, temperature, depth_of_interfaces::Number, eos)

    S_u = S_g = salinity[1]
    S_l = S_f = salinity[2]
    T_u = T_f = temperature[1]
    T_l = T_g = temperature[2]

    ρ_u = ρ(T_u, S_u, depth_of_interfaces, eos)
    ρ_l = ρ(T_l, S_l, depth_of_interfaces, eos)
    ρ_f = ρ(T_f, S_f, depth_of_interfaces, eos)
    ρ_g = ρ(T_g, S_g, depth_of_interfaces, eos)

    return (0.5 * (ρ_f - ρ_u) + 0.5 * (ρ_l - ρ_g)) / (0.5 * (ρ_f - ρ_l) + 0.5 * (ρ_u - ρ_g))
end

function Base.show(io::IO, sics::AbstractSingleInterfaceInitialConditions)
    if sics isa PeriodicSTSingleInterfaceInitialConditions
        println(io, "PeriodicSTSingleInterfaceInitialConditions")
        println(io, "┣━━ depth_of_interface: $(sics.depth_of_interface)")
        println(io, "┣━━━━━ salinity_values: $(sics.salinity_values)")
        println(io, "┣━━ temperature_values: $(sics.temperature_values)")
        println(io, "┣━ interface_smoothing: $(summary(sics.interface_smoothing))")
        println(io, "┣━━━━ background_state: $(summary(sics.background_state))")
        print(io,   "┗━━━━━━━━━━━━━━━━━ R_ρ: $(round.(sics.R_ρ; digits = 2))")
    elseif sics isa STSingleInterfaceInitialConditions
        println(io, "STSingleInterfaceInitialConditions")
        println(io, "┣━━ depth_of_interface: $(sics.depth_of_interface)")
        println(io, "┣━━━━━ salinity_values: $(sics.salinity_values)")
        println(io, "┣━━ temperature_values: $(sics.temperature_values)")
        println(io, "┣━ interface_smoothing: $(summary(sics.interface_smoothing))")
        print(io,   "┗━━━━━━━━━━━━━━━━━ R_ρ: $(round.(sics.R_ρ; digits = 2))")
    end
end
