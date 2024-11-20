abstract type AbstractStaircaseInitialConditions <: AbstractInitialConditions end

Base.iterate(sics::AbstractStaircaseInitialConditions, state = 1) =
    state > length(fieldnames(sics)) ? nothing : (getfield(sics, state), state + 1)

"""
    struct STSingleInterfaceInitialConditions
Initial conditions for a single interface (i.e. two layers with uniform `S` and `T` seperated
by a step change). The property `maintain_interface` is a `Boolean` which if set to `true` will
set [reentrant_boundary_conditions](@ref) so that the interface will be maintained (by not
letting the system run down).
"""
struct STSingleInterfaceInitialConditions{T, A, B} <: AbstractStaircaseInitialConditions
    "The depth of the interface"
    depth_of_interface :: T
    "Salinity values in each layer"
       salinity_values :: A
     "Temperature values in each layer"
    temperature_values :: A
     "Initial R_ρ at the interface"
                   R_ρ :: T
    "Boolean whether or not to set reentrant boundary condtions to approximately maintain the initial
    interface gradients"
    maintain_interface :: B
end
function STSingleInterfaceInitialConditions(model, depth_of_interface, salinity, temperature;
                                            maintain_interface = false)

    eos = model.buoyancy.model.equation_of_state

    R_ρ = compute_R_ρ(salinity, temperature, depth_of_interface, eos)

    return STSingleInterfaceInitialConditions(depth_of_interface, salinity, temperature, R_ρ, maintain_interface)

end
function STSingleInterfaceInitialConditions(eos::BoussinesqEquationOfState, depth_of_interface, salinity, temperature;
                                            maintain_interface = false)

    R_ρ = compute_R_ρ(salinity, temperature, depth_of_interface, eos)

    return STSingleInterfaceInitialConditions(depth_of_interface, salinity, temperature, R_ρ, maintain_interface)

end
const SingleInterfaceICs = STSingleInterfaceInitialConditions # alias
Base.summary(ics::SingleInterfaceICs) = "Single S-T interface at z = $(ics.depth_of_interface)"
"""
    struct PeriodicSTSingleInterfaceInitialConditions
Sets a `BackgroundField` according to `background_State` and uses a triply periodic domain
to evolve salinity and temperature anomalies about the background state.
"""
struct PeriodicSTSingleInterfaceInitialConditions{T, A, F} <: AbstractStaircaseInitialConditions
    "The depth of the interface"
    depth_of_interface :: T
    "Salinity values in each layer"
       salinity_values :: A
     "Temperature values in each layer"
    temperature_values :: A
     "Initial R_ρ at the interface"
                   R_ρ :: T
    "`BackgroundField` about which anomaly is advected. Should be an `AbstractBackgroundFunction`."
      background_state :: F
end
function PeriodicSTSingleInterfaceInitialConditions(eos::BoussinesqEquationOfState, depth_of_interface, salinity, temperature, background_state)

    R_ρ = compute_R_ρ(salinity, temperature, depth_of_interface, eos)

    return PeriodicSTSingleInterfaceInitialConditions(depth_of_interface, salinity, temperature, R_ρ, background_state)

end
const PeriodoicSingleInterfaceICs = PeriodicSTSingleInterfaceInitialConditions # alias
Base.summary(ics::PeriodoicSingleInterfaceICs) = "Single S-T interface at z = $(ics.depth_of_interface) on triply periodic domain with $(summary(ics.background_state)) state"

"Container for initial conditions that have well mixed layers seperated by sharp step interfaces."
struct STStaircaseInitialConditions{T} <: AbstractStaircaseInitialConditions
    "Number of interfaces in the initial state"
    number_of_interfaces :: Int
    "The depth of the interfaces, note length(depth_of_interfaces) == number_of_interfaces"
     depth_of_interfaces :: T
    "Salinity values in each layer"
         salinity_values :: T
    "Temperature values in each layer"
      temperature_values :: T
    "Initial R_ρ at each step interface"
                     R_ρ :: T
end
function STStaircaseInitialConditions(model, number_of_interfaces, depth_of_interfaces, salinity, temperature)

    eos = model.buoyancy.model.equation_of_state

    R_ρ = compute_R_ρ(salinity, temperature, depth_of_interfaces, eos)

    return STStaircaseInitialConditions(number_of_interfaces, depth_of_interfaces, salinity, temperature, R_ρ)

end
const StaircaseICs = STStaircaseInitialConditions # alias
Base.summary(ics::StaircaseICs) = "Multiple S-T interfaces at z = $(ics.depth_of_interfaces)"

"""
    function compute_R_ρ(salinity, temperature, depth_of_interfaces, eos)
Compute the density ratio, ``R_{\rho}``, at a diffusive interface with a non-linear equation of state
as defined in [McDougall (1981)](https://www.sciencedirect.com/science/article/abs/pii/0079661181900021)
equation (1) on page 92.
"""
function compute_R_ρ(salinity, temperature, depth_of_interfaces::Array, eos)

    R_ρ = similar(depth_of_interfaces)

    S_u = S_g = @view salinity[1:end-1]
    S_l = S_f = @view salinity[2:end]
    T_u = T_f = @view temperature[1:end-1]
    T_l = T_g = @view temperature[2:end]

    eos_vec = fill(eos, length(S_u))

    ρ_u = ρ.(T_u, S_u, depth_of_interfaces, eos_vec)
    ρ_l = ρ.(T_l, S_l, depth_of_interfaces, eos_vec)
    ρ_f = ρ.(T_f, S_f, depth_of_interfaces, eos_vec)
    ρ_g = ρ.(T_g, S_g, depth_of_interfaces, eos_vec)

    R_ρ = @. (0.5 * (ρ_f - ρ_u) + 0.5 * (ρ_l - ρ_g)) / (0.5 * (ρ_f - ρ_l) + 0.5 * (ρ_u - ρ_g))

    return R_ρ
end
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

# TODO: implement this so it is an option.
"Container for initial conditions that have well mixed layers seperated by smoothed step interfaces."
struct SmoothSTStaircaseInitialConditions{T} <: AbstractStaircaseInitialConditions
    "Number of interfaces in the initial state"
    number_of_interfaces :: Int
    "The depth of the interfaces, **note:** length (depth_of_interfaces) == number_of_interfaces"
     depth_of_interfaces :: T
    "Salinity values in each layer"
         salinity_values :: T
    "Temperature values in each layer"
      temperature_values :: T
    "Function to smooth step transition"
      smoothing_funciton :: Function
end
function Base.show(io::IO, sics::AbstractStaircaseInitialConditions)
    if sics isa PeriodicSTSingleInterfaceInitialConditions
        println(io, "STSingleInterfaceInitialConditions")
        println(io, "┣━ depth_of_interface: $(sics.depth_of_interface)")
        println(io, "┣━━━━ salinity_values: $(sics.salinity_values)")
        println(io, "┣━ temperature_values: $(sics.temperature_values)")
        println(io, "┣━━━━━━━━━━━━━━━━ R_ρ: $(round.(sics.R_ρ; digits = 2))")
        println(io, "┗━━━ background_state: $(summary(sics.background_state))")
    elseif sics isa STSingleInterfaceInitialConditions
        println(io, "STSingleInterfaceInitialConditions")
        println(io, "┣━ depth_of_interface: $(sics.depth_of_interface)")
        println(io, "┣━━━━ salinity_values: $(sics.salinity_values)")
        println(io, "┣━ temperature_values: $(sics.temperature_values)")
        println(io, "┣━ maintain_interface: $(round.(sics.maintain_interface))")
        print(io,   "┗━━━━━━━━━━━━━━━━ R_ρ: $(round.(sics.R_ρ; digits = 2))")
    elseif sics isa STStaircaseInitialConditions
        println(io, "STStaircaseInitialConditions")
        println(io, "┣━ number_of_interfaces: $(sics.number_of_interfaces)")
        println(io, "┣━━ depth_of_interfaces: $(sics.depth_of_interfaces)")
        println(io, "┣━━━━━━ salinity_values: $(sics.salinity_values)")
        println(io, "┣━━━ temperature_values: $(sics.temperature_values)")
        print(io,   "┗━━━━━━━━━━━━━━━━━━ R_ρ: $(round.(sics.R_ρ; digits = 2))")
    elseif sics isa SmoothSTStaircaseInitialConditions
        println(io, "STStaircaseInitialConditions")
        println(io, "┣━ number_of_interfaces: $(sics.number_of_interfaces)")
        println(io, "┣━━ depth_of_interfaces: $(sics.depth_of_interfaces)")
        println(io, "┣━━━━━━ salinity_values: $(sics.salinity_values)")
        println(io, "┣━━━ temperature_values: $(sics.temperature_values)")
        print(io,   "┗━━━ smoothing_funciton: $(sics.smoothing_funciton)")
    end
end
