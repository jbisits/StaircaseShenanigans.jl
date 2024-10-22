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
@inline (quarter::OuterMask)(x, y, z) = z > quarter.upper ? 1 : z < quarter.lower ? 1 : 0

Base.summary(oqm::OuterMask) = "Upper masked region [$(oqm.upper), 0), lower masked region [-Lz, $(oqm.lower))"

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
Base.summary(oqt::OuterMask) = "Upper target = $(oqt.upper), lower target = $(oqt.lower)"

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
