# The following functions are to be used as `BoundaryConditions` so that tracers can
# re-enter the domain with the initial gradient added effectively allowing the gradient to
# be maintained. This is effectively a restoring (or more maintaining) which is why it
# is in this script.

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
