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
OuterStairMask(depth_of_steps) = OuterStairMask(depth_of_steps[1], depth_of_steps[end])

@inline (stairs::OuterStairMask)(x, y, z) =
            z > stairs.first_stair ? 1 : z < stairs.last_stair ? 1 : 0

Base.summary(oms::OuterStairMask) = "First stair mask [0, $(oms.first_stair)], last stair mask [$(oms.last_stair), end]"
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
OuterStairTargets(depth_of_steps, initial_tracer_values) =
    OuterStairTargets(depth_of_steps[1], initial_tracer_values[1],
                      depth_of_steps[end], initial_tracer_values[end])

@inline (stair_targets::OuterStairTargets)(x, y, z, t) =
            z > stair_targets.first_stair ? stair_targets.first_target :
            z < stair_targets.last_stair ? stair_targets.last_target : 0

Base.summary(ost::OuterStairTargets) = "First stair target $(ost.first_target), last stair target $(ost.last_target)"
