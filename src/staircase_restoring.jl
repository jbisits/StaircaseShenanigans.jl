Oceananigans.BackgroundField(bf::BackgroundFunction) =
    BackgroundField(bf.func, parameters = bf.parameters)
"""
    mutable struct BackgroundTanh{F, T, P}
Container for a tanh background field.
"""
mutable struct BackgroundTanh{F, T} <: BackgroundFunction
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
mutable struct BackgroundLinear{F} <: BackgroundFunction
    "Linear function"
          func :: F
    "Parameters for the tanh background field"
    parameters :: NamedTuple
end
BackgroundLinear() = BackgroundLinear(linear_background, NamedTuple())

function Base.show(io::IO, bf::BackgroundFunction)
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
Base.summary(bt::BackgroundTanh) = "$(bt.func)"
Base.summary(bl::BackgroundLinear) = "$(bl.func)"

"""
    function S_and_T_background_fields(initial_conditions)
Set background fields for the `S` and `T` tracer fields where the domain is triply periodic.
"""
function S_and_T_background_fields(ics::PeriodicSTSingleInterfaceInitialConditions, Lz)

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

function get_parameters!(ics::PeriodicSTSingleInterfaceInitialConditions, tracer::Symbol, Lz)

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
function save_background_state!(simulation, model, initial_conditions::PeriodoicSingleInterfaceICs)

    S_background = Field(model.background_fields.tracers.S)
    compute!(S_background)
    T_background = Field(model.background_fields.tracers.T)
    compute!(T_background)
    σ_background = Field(seawater_density(model, temperature = T_background, salinity = S_background,
                                         geopotential_height = 0))
    compute!(σ_background)

    if simulation.output_writers[:tracers] isa NetCDFOutputWriter

        NCDataset(simulation.output_writers[:tracers].filepath, "a") do ds
            defVar(ds, "S_background", S_background, ("xC", "yC", "zC"),
                  attrib = Dict("longname" => "Background field for salinity",
                                "units" => "gkg⁻¹"))
            defVar(ds, "T_background", T_background, ("xC", "yC", "zC"),
                  attrib =  Dict("longname" => "Background field for temperature",
                                 "units" => "°C"))
        end

        NCDataset(simulation.output_writers[:computed_output].filepath, "a") do ds
            defVar(ds, "σ_background", σ_background, ("xC", "yC", "zC"),
                  attrib = Dict("longname" => "Background field for potential density (0dbar) computed from the `S` and `T` background fields",
                                "units" => "kgm⁻³"))
        end

    elseif simulation.output_writers[:tracers] isa JLD2OutputWriter

        jldopen(simulation.output_writers[:tracers].filepath, "a+") do f
            f["S_background"] = S_background
            f["T_background"] = T_background
        end

        jldopen(simulation.output_writers[:computed_output].filepath, "a+") do f
            f["σ_background"] = σ_background
        end

    end

    return nothing
end

# The following functions are to be used as `BoundaryConditions` so that tracers can
# re-enter the domain with the initial gradient added effectively allowing the gradient to
# be maintained. Another option is to add background tracer fields and only evolve the anomaly

"""
    T_reentry(i, j, grid, clock, model_fields, ΔT)
Discrete boundary condition to set the temperature on a vertical boundary to tracer value at
the **top** of the domain and adds ΔT. That is a reentrant T condition to be used on the
**bottom** of the domain.
"""
@inline T_reentry(i, j, grid, clock, model_fields, ΔT) =
    @inbounds ℑzᵃᵃᶠ(i, j, grid.Nz+1, model_fields.T) + ΔT
"""
    S_reentry(i, j, grid, clock, model_fields, ΔS)
As for [T_reentry](@ref) but using salinity tracer instead.
"""
@inline S_reentry(i, j, grid, clock, model_fields, ΔS) =
    @inbounds ℑzᵃᵃᶠ(i, j, grid.Nz+1, model_fields.S) + ΔS

"Access velocity field at the _bottom_ of the domain for use as `BoundaryCondition`."
@inline _w_bottom(i, j, grid, clock, model_fields) = @inbounds model_fields.w[i, j, 1]
"Access velocity field at the _top_ of the domain for use as `BoundaryCondition`."
@inline _w_top(i, j, grid, clock, model_fields) = @inbounds model_fields.w[i, j, grid.Nz+1]

"Interpolated velocity between top and bottom (vertical) faces of domain."
@inline w_top_bottom_interpolate(i, j, grid, clock, model_fields) =
    @inbounds 0.5 * (model_fields.w[i, j, 1] + model_fields.w[i, j, grid.Nz+1])

""""
    function reentrant_boundary_conditions(ics::SingleInterfaceICs)
Setup boundary conditions to maintain a gradient (ΔC) between the two layers of a
`SingleInterface` model. The boundary conditions for the tracers are re-entrant from top to
bottom with ΔC added to maintain a ΔC difference between the top and bottom layer. The
vertical velocitiy field also has modified vertical boundary conditions where the velocity
on the bottom face is set to the velocity on the top face and vice versa.
"""
function reentrant_boundary_conditions(ics::SingleInterfaceICs)

    bcs = if ics.maintain_interface
            ΔT = diff(ics.temperature_values)[1]
            ΔS = diff(ics.salinity_values)[1]
            T_bottom_reentry = ValueBoundaryCondition(T_reentry, discrete_form=true, parameters = ΔT)
            S_bottom_reentry = ValueBoundaryCondition(S_reentry, discrete_form=true, parameters = ΔS)
            T_bcs = FieldBoundaryConditions(bottom = T_bottom_reentry)
            S_bcs = FieldBoundaryConditions(bottom = S_bottom_reentry)

            w_top = OpenBoundaryCondition(_w_bottom, discrete_form=true)
            w_bottom = OpenBoundaryCondition(_w_top, discrete_form=true)
            w_bcs = FieldBoundaryConditions(top = w_top, bottom = w_bottom)

            (T=T_bcs, S=S_bcs, w=w_bcs)
        else
            NamedTuple()
        end

    return bcs
end
"""
    struct OuterStairMask{T}
Container for location of first and last stair for `mask` in `Relaxation`. This need not be
the whole mixed region above/below the first/last stair. It can be
"""
struct OuterStairMask{T}
    "location of the first stair in the thermohaline staircase"
    first_stair :: T
    "location of the last stair in the thermohaline staircase"
     last_stair :: T
end
"Convenience constructor if depth of steps is defined for use in initial condition setting"
OuterStairMask(depth_of_interfaces) = OuterStairMask(depth_of_interfaces[1], depth_of_interfaces[end])

@inline (stairs::OuterStairMask)(x, y, z) =
            z > stairs.first_stair ? 1 : z < stairs.last_stair ? 1 : 0

Base.summary(oms::OuterStairMask) = "First stair mask [$(oms.first_stair), 0], last stair mask [-Lz, $(oms.last_stair)]"
"""
    struct OuterStairTargets{T}
Location and target values for `target` use in `Relaxation`.
"""
struct OuterStairTargets{T}
    "location of the first stair in the thermohaline staircase"
    first_stair :: T
    "target value for the first stair in the thermohaline staircase"
    first_target :: T
    "location of the last stair in the thermohaline staircase"
    last_stair :: T
    "target value for the first stair in the thermohaline staircase"
    last_target :: T
end
"Convenience constructor for restoring to initial tracer values"
OuterStairTargets(depth_of_interfaces, initial_tracer_values) =
    OuterStairTargets(depth_of_interfaces[1], initial_tracer_values[1],
                      depth_of_interfaces[end], initial_tracer_values[end])

@inline (stair_targets::OuterStairTargets)(x, y, z, t) =
            z > stair_targets.first_stair ? stair_targets.first_target :
            z < stair_targets.last_stair ? stair_targets.last_target : 0

Base.summary(ost::OuterStairTargets) = "First stair target $(ost.first_target), last stair target $(ost.last_target)"

"""
    struct OuterMask{T}
Container with the upper and lower regions of of the domain to mask.
"""
struct OuterMask{T}
    "End of upper quarter of domain, z ∈ [-upper, 0)"
    upper :: T
    "Start of lower quarter of domain, z ∈ [-Lz, lower)"
    lower :: T
end

"Function to find upper and lower regions of domain."
@inline (om::OuterMask)(x, y, z) = z > om.upper ? 1 : z < om.lower ? 1 : 0

Base.summary(om::OuterMask) = "Upper masked region [$(om.upper), 0), lower masked region [-Lz, $(om.lower))"

"""
    struct OuterTargets{T}
Target values for the upper and lower regions of domain
"""
struct OuterTargets{T}
    "Upper region target value"
    upper_target :: T
    "Lower region target value"
    lower_target :: T
    "Upper region depth"
    upper_depth :: T
end
Base.summary(oqt::OuterTargets) = "Upper target = $(oqt.upper_target), lower target = $(oqt.lower_target)"

"Convenience constructor to get the upper region of domain from an `OuterMask`"
OuterTargets(upper_target, lower_target, mask::OuterMask) =
    OuterTargets(upper_target, lower_target, mask.upper)
"Restoring function"
@inline (oqt::OuterTargets)(x, y, z, t) = z > oqt.upper_depth ? oqt.upper_target :
                                                                oqt.lower_target

"""
    struct ExponentialTarget{T}
Restore to the tracer ``C`` to exponetial function of the form ``C(x, y, z, t) = Aℯ⁻ᶻ``.
"""
struct ExponentialTarget{T}
    A :: T
    λ :: T
end

@inline (p::ExponentialTarget)(x, y, z, t) = p.A * exp(-p.λ * z)
