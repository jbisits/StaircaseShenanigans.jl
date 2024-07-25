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
function StepInitialConditions(model, number_of_steps, depth_of_steps, salinity, temperature)

    grid = model.grid
    eos = model.buoyancy.model.equation_of_state
    ΔT = diff(temperature)
    ΔS = diff(salinity)
    α, β = compute_α_and_β(salinity, temperature, depth_of_steps, eos, grid)

    R_ρ = similar(depth_of_steps)
    @. R_ρ = (β * ΔS) / (α * ΔT)

    return StepInitialConditions(number_of_steps, depth_of_steps, salinity, temperature, R_ρ)

end
"""
    function compute_α_and_β(salinity, temperature, eos)
Compute thermal expansion and haline contraction coefficients at the interface of the steps.
The coefficients are computed as α̂ = 0.5 * (α(Sᵢ, 0.5 * (Θᵢ+Θⱼ), 0) + α(Sⱼ, 0.5 * (Θᵢ+Θⱼ), 0))
where j = i + 1. Still need to double check if there is a more accurate way to do this as
I think this is a slight simplication of the method I am meant to be using.
"""
function compute_α_and_β(salinity, temperature, depth_of_steps, eos, grid)

    S1 = @view salinity[1:end-1]
    S2 = @view salinity[2:end]
    Sstepmean = (S1 .+ S2) / 2

    T1 = @view temperature[1:end-1]
    T2 = @view temperature[2:end]
    Tstepmean = (T1 .+ T2) / 2

    eos_vec = fill(eos, length(S1))

    gh = Field(geopotential_height(grid))
    compute!(gh)

    step_gh_idx = [findfirst(interior(gh, 1, 1, :) .≥ d) for d ∈ depth_of_steps]
    step_gh_height = interior(gh, 1, 1, step_gh_idx)

    α = 0.5 * (thermal_expansion.(Tstepmean, S1, step_gh_height, eos_vec) +
               thermal_expansion.(Tstepmean, S2, step_gh_height, eos_vec))

    β = 0.5 * (haline_contraction.(T1, Sstepmean, step_gh_height, eos_vec) +
               haline_contraction.(T2, Sstepmean, step_gh_height, eos_vec))

    return α, β
end
geopotential_height(grid) = KernelFunctionOperation{Center, Center, Face}(Zᶜᶜᶠ, grid)

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
