Oceananigans.BackgroundField(bf::AbstractBackgroundFunction) =
    BackgroundField(bf.func, parameters = bf.parameters)

"""
    mutable struct BackgroundTanh{F, T, P}
Container for a tanh background field.
"""
mutable struct BackgroundTanh{F, T} <: AbstractBackgroundFunction
    "tanh function"
          func :: F
    "Thickness of salinity interface"
            hₛ :: T
    "Thickness of temperature interface"
            hₜ :: T
    "Parameters for the tanh background field"
    parameters :: NamedTuple
end
BackgroundTanh() = BackgroundTanh(tanh_background, 0.02, 0.04, NamedTuple())
"Default that makes the temperature interface 3 times larger than the salinity"
function BackgroundTanh(scale::Number)

    hₛ = scale
    hₜ = 3 * scale
    _hₛ, _hₜ = promote(hₛ, hₜ) # ensures both are same type

    return BackgroundTanh(tanh_background, _hₛ, _hₜ, NamedTuple())
end
"Pass a `NamedTuple` in form `(S = , T = )` to set the tanh scaling."
function BackgroundTanh(scale::NamedTuple)

    hₛ, hₜ = scale.S, scale.T
    _hₛ, _hₜ = promote(hₛ, hₜ) # ensures both are same type

    return BackgroundTanh(tanh_background, _hₛ, _hₜ, NamedTuple())
end
"""
    mutable struct BackgroundLinear{F}
Container for a linear background field.
"""
mutable struct BackgroundLinear{F} <: AbstractBackgroundFunction
    "Linear function"
          func :: F
    "Parameters for the linear background field"
    parameters :: NamedTuple
end
BackgroundLinear() = BackgroundLinear(linear_background, NamedTuple())

"""
    mutable struct BackgroundStep{F}
Container for a step change background field.
"""
mutable struct BackgroundStep{F} <: AbstractBackgroundFunction
    "Heaviside function"
          func :: F
    "Parameters for the step change background field"
    parameters :: NamedTuple
end
BackgroundStep() = BackgroundStep(step_background, NamedTuple())

function Base.show(io::IO, bf::AbstractBackgroundFunction)
    if bf isa BackgroundTanh
        println(io, "BackgroundTanh")
        println(io, "┣━━━ function: $(bf.func)")
        println(io, "┣━━━━ S scale: $(bf.hₛ)")
        println(io, "┣━━━━ T scale: $(bf.hₜ)")
          print(io, "┗━ parameters: $(bf.parameters)")
    else
        println(io, "$(typeof(bf))")
        println(io, "┣━━━ function: $(bf.func)")
          print(io, "┗━ parameters: $(bf.parameters)")
    end
end

struct NoBackgroundFunction <: AbstractBackgroundFunction end
const NoBackground = NoBackgroundFunction
Base.summary(bt::BackgroundTanh) = "$(bt.func)"
Base.summary(bl::BackgroundLinear) = "$(bl.func)"
Base.summary(bn::Type{<:NoBackgroundFunction}) = "no background field"

S_and_T_background_fields(ics::STSingleInterfaceInitialConditions, z_range) =
    S_and_T_background_fields(ics, z_range, ics.background_state)
"Return blank `NamedTuple` so that no `BackgroundField`s are set."
S_and_T_background_fields(ics, z_range, background_state::Type{<:NoBackground}) = NamedTuple()
"""
    function S_and_T_background_fields(initial_conditions)
Set background fields for the `S` and `T` tracer fields where the domain is triply periodic.
"""
function S_and_T_background_fields(ics, z_range, background_state)

    get_parameters!(ics, :salinity_values, z_range)
    S_background = BackgroundField(ics.background_state)
    get_parameters!(ics, :temperature_values, z_range)
    T_background = BackgroundField(ics.background_state)

    return (S = S_background, T = T_background)
end

S_and_T_background_fields(ics, z_range, background_state::NamedTuple) = background_state

"Sets a background state that is hyperbolic tangent. There is also a method to save an
`Array` of this backgorund state to output."
@inline tanh_background(x, y, z, t, p) = p.Cₗ - 0.5 * p.ΔC * (1  + tanh( (z - p.z_interface) / (p.h * p.Δz)))
"Sets a background state that is linear. There is also a method to save an
`Array` of this backgorund state to output."
@inline linear_background(x, y, z, t, p) = p.Cᵤ + (p.ΔC / p.Δz) * (z - p.z_top)
"Sets a background state that is step function. There is also a method to save an
`Array` of this backgorund state to output."
@inline step_background(x, y, z, t, p) = z > p.z_interface ? p.Cᵤ : p.Cₗ

function get_parameters!(ics::STSingleInterfaceInitialConditions, tracer::Symbol, z_range)

    z_interface = ics.depth_of_interface
    C = Array(getproperty(ics, tracer))
    ΔC = diff(C)[1]
    Cᵤ, Cₗ = C
    Lz, z_top =  z_range
    Δz = Lz - z_top
    update_parameters!(ics.background_state, tracer, ΔC, Cᵤ, Cₗ, Δz, z_top, z_interface)

    return nothing
end
function update_parameters!(background_state::BackgroundTanh, tracer, ΔC, Cᵤ, Cₗ, Δz, z_top, z_interface)

    h = tracer == :salinity_values ? background_state.hₛ : background_state.hₜ
    background_state.parameters = (; ΔC, Cᵤ, Cₗ, Δz = abs(Δz), z_interface, h)

    return nothing
end
function update_parameters!(background_state::BackgroundLinear, tracer, ΔC, Cᵤ, Cₗ, Δz, z_top, z_interface)

    background_state.parameters = (; ΔC, Δz, Cᵤ, z_top)

    return nothing
end
function update_parameters!(background_state::BackgroundStep, tracer, ΔC, Cᵤ, Cₗ, Δz, z_top, z_interface)

    background_state.parameters = (; Cᵤ, Cₗ, z_interface)

    return nothing
end

"""
    function save_background_state!(simulation, sdns)
Where there is `BackgroundField` (currently this assumes that is for periodic simulations)
save the background state for the tracers and density so this can be used later
"""
save_background_state!(simulation, sdns) = save_background_state!(simulation, sdns.model, sdns.initial_conditions.background_state)
save_background_state!(simulation, model, background_state::Type{<:NoBackground}) = nothing
function save_background_state!(simulation, model, background_state)

    S_background = Field(model.background_fields.tracers.S)
    compute!(S_background)
    S_background_array = Array(interior(S_background, :, :, :))
    T_background = Field(model.background_fields.tracers.T)
    compute!(T_background)
    T_background_array = Array(interior(T_background, :, :, :))
    σ_background = Field(seawater_density(model, temperature = T_background, salinity = S_background,
                                        geopotential_height = 0))
    compute!(σ_background)
    σ_background_array = Array(interior(σ_background, :, :, :))

    # Here need only check if one key is present as if one is all are
    if !check_key_present(simulation, :tracers, "S_background")

        ow = simulation.output_writers
        if ow[:tracers] isa NetCDFWriter

            NCDataset(ow[:tracers].filepath, "a") do ds
                defVar(ds, "S_background", S_background_array, ("xC", "yC", "zC"),
                    attrib = Dict("longname" => "Background field for salinity",
                                    "units" => "gkg⁻¹"))
                defVar(ds, "T_background", T_background_array, ("xC", "yC", "zC"),
                    attrib =  Dict("longname" => "Background field for temperature",
                                    "units" => "°C"))
            end

            NCDataset(ow[:computed_output].filepath, "a") do ds
                defVar(ds, "σ_background", σ_background_array, ("xC", "yC", "zC"),
                    attrib = Dict("longname" => "Background field for potential density (0dbar) computed from the `S` and `T` background fields",
                                    "units" => "kgm⁻³"))
            end

        elseif ow[:tracers] isa JLD2Writer

            jldopen(ow[:tracers].filepath, "a+") do f
                f["S_background"] = S_background_array
                f["T_background"] = T_background_array
            end

            jldopen(ow[:computed_output].filepath, "a+") do f
                f["σ_background"] = σ_background_array
            end

        end

    end
    return nothing
end
