abstract type AbstractStaircaseInitialConditions <: AbstractInitialConditions end

Base.iterate(sics::AbstractStaircaseInitialConditions, state = 1) =
    state > length(fieldnames(sics)) ? nothing : (getfield(sics, state), state + 1)

"Container for initial conditions that have well mixed layers seperated by sharp step interfaces."
struct STStaircaseInitialConditions{I, T, IS, BF} <: AbstractStaircaseInitialConditions
    "Number of interfaces in the initial state"
    number_of_interfaces :: I
    "The depth of the interfaces, note length(depth_of_interfaces) == number_of_interfaces"
     depth_of_interfaces :: T
    "Salinity values in each layer"
         salinity_values :: T
    "Temperature values in each layer"
      temperature_values :: T
      "Initial smoothing function overeach interface"
     interface_smoothing :: IS
      "`BackgroundField` about which anomaly is advected. Should be an `AbstractBackgroundFunction`."
        background_state :: BF
    "Initial R_ρ at each step interface"
                     R_ρ :: T
end
function STStaircaseInitialConditions(model, number_of_interfaces, depth_of_interfaces,
                                      salinity, temperature;
                                      interface_smoothing = NoSmoothing,
                                      background_state = NoBackground)

    eos = model.buoyancy.formulation.equation_of_state

    R_ρ = compute_R_ρ(salinity, temperature, depth_of_interfaces, eos)

    return STStaircaseInitialConditions(number_of_interfaces, depth_of_interfaces, salinity,
                                        temperature, interface_smoothing, background_state, R_ρ)

end
function STStaircaseInitialConditions(eos::BoussinesqEquationOfState, number_of_interfaces,
                                      depth_of_interfaces, salinity, temperature;
                                      interface_smoothing = NoSmoothing,
                                      background_state = NoBackground)

    R_ρ = compute_R_ρ(salinity, temperature, depth_of_interfaces, eos)

    return STStaircaseInitialConditions(number_of_interfaces, depth_of_interfaces, salinity,
                                        temperature, interface_smoothing, background_state, R_ρ)

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

    ρ_u = total_density.(T_u, S_u, depth_of_interfaces, eos_vec)
    ρ_l = total_density.(T_l, S_l, depth_of_interfaces, eos_vec)
    ρ_f = total_density.(T_f, S_f, depth_of_interfaces, eos_vec)
    ρ_g = total_density.(T_g, S_g, depth_of_interfaces, eos_vec)

    R_ρ = @. (0.5 * (ρ_f - ρ_u) + 0.5 * (ρ_l - ρ_g)) / (0.5 * (ρ_f - ρ_l) + 0.5 * (ρ_u - ρ_g))

    return R_ρ
end

function Base.show(io::IO, sics::AbstractStaircaseInitialConditions)
    println(io, "STStaircaseInitialConditions")
    println(io, "┣━ number_of_interfaces: $(sics.number_of_interfaces)")
    println(io, "┣━━ depth_of_interfaces: $(sics.depth_of_interfaces)")
    println(io, "┣━━━━━━ salinity_values: $(sics.salinity_values)")
    println(io, "┣━━━ temperature_values: $(sics.temperature_values)")
    println(io, "┣━ interface_smoothing: $(summary(sics.interface_smoothing))")
    println(io, "┣━━━━ background_state: $(summary(sics.background_state))")
    print(io,   "┗━━━━━━━━━━━━━━━━━━ R_ρ: $(round.(sics.R_ρ; digits = 2))")
end
