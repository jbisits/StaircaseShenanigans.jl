Base.iterate(sics::AbstractStaircaseInitialConditions, state = 1) =
    state > length(fieldnames(sics)) ? nothing : (getfield(sdns, state), state + 1)

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
function StepInitialConditions(number_of_steps, depth_of_steps, salinity, temperature)

    ΔT = diff(temperature)
    ΔS = diff(salinity)
    # Constant values take from Vallis, not sure if this is the best way to calculate
    # these but will start with this.
    α = 1.67e-4
    β = 7.80e-4

    R_ρ = similar(depth_of_steps)
    @. R_ρ = (β * ΔS) / (α * ΔT)

    return StepInitialConditions(number_of_steps, depth_of_steps, salinity, temperature, R_ρ)

end
"""
    function set_steps(z, C, depth_of_steps)
Set the salinity and temperature step initial conditions. The salinity and temperature are
set in the ranges:
- z ∈ [depth_of_steps[1], 0)
- z ∈ [depth_of_steps[i], depth_of_steps[i-1]) for i = 2:number_of_steps-1
- z ∈ [-Lz, depth_of_steps[1]).

**Not currently in use and will remove if not needed.**
"""
function set_steps(z, C, depth_of_steps)

    # This might need to be set as an `Array` rather than a function as at the moment
    # I cannot think of a way to generalise this (I have tried using loops but it breaks..).
    # For now I am goint to experiment with some four step experiments and see how things go.
    if z > depth_of_steps[1]
        C[1]
    elseif depth_of_steps[2] < z ≤ depth_of_steps[1]
        C[2]
    elseif depth_of_steps[3] < z ≤ depth_of_steps[2]
        C[3]
    elseif depth_of_steps[4] < z ≤ depth_of_steps[3]
        C[4]
    elseif z < depth_of_steps[4]
        C[5]
    end

end
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
