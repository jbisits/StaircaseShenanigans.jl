abstract type AbstractStaircaseInitialConditions <: AbstractInitialConditions end

Base.iterate(sics::AbstractStaircaseInitialConditions, state = 1) =
    state > length(fieldnames(sics)) ? nothing : (getfield(sics, state), state + 1)

"Container for initial conditions that have well mixed layers seperated by sharp step interfaces."
struct StepInitialConditions{T} <: AbstractStaircaseInitialConditions
    "Number of step changes in the initial state"
       number_of_steps :: Int
    "The depth at which the step changes take place, note length(depth_of_steps) == number_of_steps"
        depth_of_steps :: T
    "Salinity values in each layer"
       salinity_values :: T
    "Temperature values in each layer"
    temperature_values :: T
    "Initial R_ρ at each step interface"
                   R_ρ :: T
end
function StepInitialConditions(model, number_of_steps, depth_of_steps, salinity, temperature)

    eos = model.buoyancy.model.equation_of_state

    R_ρ = compute_R_ρ(salinity, temperature, depth_of_steps, eos)

    return StepInitialConditions(number_of_steps, depth_of_steps, salinity, temperature, R_ρ)

end
"""
    function compute_R_ρ(salinity, temperature, depth_of_steps, eos)
Compute the density ratio, ``R_{\rho}``, at a diffusive interface with a non-linear equation of state
as defined in [McDougall (1981)](https://www.sciencedirect.com/science/article/abs/pii/0079661181900021)
equation (1) on page 92.
"""
function compute_R_ρ(salinity, temperature, depth_of_steps, eos)

    R_ρ = similar(depth_of_steps)

    S_u = S_g = @view salinity[1:end-1]
    S_l = S_f = @view salinity[2:end]
    Θ_u = Θ_f = @view temperature[1:end-1]
    Θ_l = Θ_g = @view temperature[2:end]

    eos_vec = fill(eos, length(S_u))

    ρ_u = ρ.(Θ_u, S_u, depth_of_steps, eos_vec)
    ρ_l = ρ.(Θ_l, S_l, depth_of_steps, eos_vec)
    ρ_f = ρ.(Θ_f, S_f, depth_of_steps, eos_vec)
    ρ_g = ρ.(Θ_g, S_g, depth_of_steps, eos_vec)

    R_ρ = @. (0.5 * (ρ_f - ρ_u) + 0.5 * (ρ_l - ρ_g)) / (0.5 * (ρ_f - ρ_l) + 0.5 * (ρ_u - ρ_g))

    return R_ρ
end

"Container for initial conditions that have well mixed layers seperated by smoothed step interfaces."
struct SmoothStepInitialConditions{T} <: AbstractStaircaseInitialConditions
    "Number of step changes in the initial state"
       number_of_steps :: Int
    "The depth at which the step changes take place, **note:** length (depth_of_steps) == number_of_steps"
        depth_of_steps :: T
    "Salinity values in each layer"
       salinity_values :: T
    "Temperature values in each layer"
    temperature_values :: T
    "Function to smooth step transition"
    smoothing_funciton :: Function
end
function Base.show(io::IO, sics::AbstractStaircaseInitialConditions)
    if sics isa StepInitialConditions
        println(io, "StepInitialConditions")
        println(io, "┣━━━━ number_of_steps: $(sics.number_of_steps)")
        println(io, "┣━━━━━ depth_of_steps: $(sics.depth_of_steps)")
        println(io, "┣━━━━ salinity_values: $(sics.salinity_values)")
        println(io, "┣━ temperature_values: $(sics.temperature_values)")
        print(io,   "┗━━━━━━━━━━━━━━━━ R_ρ: $(round.(sics.R_ρ; digits = 2))")
    elseif sics isa SmoothStepInitialConditions
        println(io, "StepInitialConditions")
        println(io, "┣━━━━ number_of_steps: $(sics.number_of_steps)")
        println(io, "┣━━━━━ depth_of_steps: $(sics.depth_of_steps)")
        println(io, "┣━━━━ salinity_values: $(sics.salinity_values)")
        println(io, "┣━ temperature_values: $(sics.temperature_values)")
        print(io,   "┗━ smoothing_funciton: $(sics.smoothing_funciton)")
    end
end
