Oceananigans.BackgroundField(bf::AbstractBackgroundFunction) =
    BackgroundField(bf.func, parameters = bf.parameters)
"""
    mutable struct BackgroundTanh{F, T, P}
Container for a tanh background field.
"""
mutable struct BackgroundTanh{F, T} <: AbstractBackgroundFunction
    "tanh function"
          func :: F
    "Scale the steepness of the `tanh` change"
         scale :: T
    "Parameters for the tanh background field"
    parameters :: NamedTuple
end
BackgroundTanh() = BackgroundTanh(tanh_background, 100, NamedTuple())
BackgroundTanh(scale) = BackgroundTanh(tanh_background, scale, NamedTuple())
"""
    mutable struct BackgroundLinear{F, P}
Container for a linear background field.
"""
mutable struct BackgroundLinear{F} <: AbstractBackgroundFunction
    "Linear function"
          func :: F
    "Parameters for the tanh background field"
    parameters :: NamedTuple
end
BackgroundLinear() = BackgroundLinear(linear_background, NamedTuple())

function Base.show(io::IO, bf::AbstractBackgroundFunction)
    if bf isa BackgroundTanh
        println(io, "BackgroundTanh")
        println(io, "┣━━━ function: $(bf.func)")
        println(io, "┣━━━━━━ scale: $(bf.scale)")
          print(io, "┗━ parameters: $(bf.parameters)")
    elseif bf isa BackgroundLinear
        println(io, "BackgroundLinear")
        println(io, "┣━━━ function: $(bf.func)")
          print(io, "┗━ parameters: $(bf.parameters)")
    end
end

struct NoBackgroundFunction <: AbstractBackgroundFunction end
const NoBackground = NoBackgroundFunction
Base.summary(bt::BackgroundTanh) = "$(bt.func)"
Base.summary(bl::BackgroundLinear) = "$(bl.func)"
Base.summary(bn::Type{<:NoBackgroundFunction}) = "no background field"

S_and_T_background_fields(ics::STSingleInterfaceInitialConditions, Lz) =
    S_and_T_background_fields(ics, Lz, ics.background_state)
"Return blank `NamedTuple` so that no `BackgroundField`s are set."
S_and_T_background_fields(ics, Lz, background_state::Type{<:NoBackground}) = NamedTuple()
"""
    function S_and_T_background_fields(initial_conditions)
Set background fields for the `S` and `T` tracer fields where the domain is triply periodic.
"""
function S_and_T_background_fields(ics, Lz, background_state)

    get_parameters!(ics, :salinity_values, Lz)
    S_background = BackgroundField(ics.background_state)
    get_parameters!(ics, :temperature_values, Lz)
    T_background = BackgroundField(ics.background_state)

    return (S = S_background, T = T_background)
end
"Sets a background state that is hyperbolic tangent. There is also a method to save an
`Array` of this backgorund state to output."
@inline tanh_background(x, y, z, t, p) = p.Cₗ - 0.5 * p.ΔC * (1  + tanh(p.D * (z - p.z_interface) / p.Lz))
tanh_background(z, ΔC, Cₗ, Lz, z_interface, D) = Cₗ - 0.5 * ΔC * (1  + tanh(D * (z - z_interface) / Lz))
@inline linear_background(x, y, z, t, p) = p.Cᵤ - p.ΔC * z / p.Lz
linear_background(z, ΔC, Cᵤ, Lz) = Cᵤ - ΔC * z / Lz

function get_parameters!(ics::STSingleInterfaceInitialConditions, tracer::Symbol, Lz)

    z_interface = ics.depth_of_interface
    C = Array(getproperty(ics, tracer))
    ΔC = diff(C)[1]
    Cᵤ, Cₗ = C
    update_parameters!(ics.background_state, ΔC, Cᵤ, Cₗ, abs(Lz), z_interface)

    return nothing
end

function update_parameters!(backgound_state::BackgroundTanh, ΔC, Cᵤ, Cₗ, Lz, z_interface)

    backgound_state.parameters = (; ΔC, Cᵤ, Cₗ, Lz, z_interface, D = backgound_state.scale)
    return nothing
end
function update_parameters!(backgound_state::BackgroundLinear, ΔC, Cᵤ, Cₗ, Lz, z_interface)

    backgound_state.parameters = (; ΔC, Cᵤ, Lz)
    return nothing
end

"""
    function save_background_state!(simulation, sdns)
Where there is `BackgroundField` (currently this assumes that is for periodic simulations)
save the background state for the tracers and density so this can be used later
"""
save_background_state!(simulation, sdns) = save_background_state!(simulation, sdns.model, sdns.initial_conditions)
save_background_state!(simulation, model, initial_conditions) = nothing
function save_background_state!(simulation, model, initial_conditions::SingleInterfaceICs)

    if !isnothing(initial_conditions.background_state)
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

        if simulation.output_writers[:tracers] isa NetCDFOutputWriter

            NCDataset(simulation.output_writers[:tracers].filepath, "a") do ds
                defVar(ds, "S_background", S_background_array, ("xC", "yC", "zC"),
                    attrib = Dict("longname" => "Background field for salinity",
                                    "units" => "gkg⁻¹"))
                defVar(ds, "T_background", T_background_array, ("xC", "yC", "zC"),
                    attrib =  Dict("longname" => "Background field for temperature",
                                    "units" => "°C"))
            end

            NCDataset(simulation.output_writers[:computed_output].filepath, "a") do ds
                defVar(ds, "σ_background", σ_background_array, ("xC", "yC", "zC"),
                    attrib = Dict("longname" => "Background field for potential density (0dbar) computed from the `S` and `T` background fields",
                                    "units" => "kgm⁻³"))
            end

        elseif simulation.output_writers[:tracers] isa JLD2OutputWriter

            jldopen(simulation.output_writers[:tracers].filepath, "a+") do f
                f["S_background"] = S_background_array
                f["T_background"] = T_background_array
            end

            jldopen(simulation.output_writers[:computed_output].filepath, "a+") do f
                f["σ_background"] = σ_background_array
            end

        end

    end
    return nothing
end
