"""
    function set_staircase_initial_conditions!(model, ics)
Set initial salinity and temperature staircases in `model`.
"""
function set_staircase_initial_conditions!(model, ics::StepInitialConditions)

    depth_of_steps = ics.depth_of_steps
    z = znodes(model.grid, Center())
    S = ics.salinity_values
    T = ics.temperature_values

    S₀ = similar(interior(model.tracers.S, :, :, :))
    T₀ = similar(interior(model.tracers.T, :, :, :))
    S₀[:, :, z .> depth_of_steps[1]] .= S[1]
    T₀[:, :, z .> depth_of_steps[1]] .= T[1]
    for i ∈ 1:length(depth_of_steps)-1
        S₀[:, :, depth_of_steps[i+1] .≤ z .< depth_of_steps[i]] .= S[i]
        T₀[:, :, depth_of_steps[i+1] .≤ z .< depth_of_steps[i]] .= T[i]
    end
    S₀[:, :, z .< depth_of_steps[end]] .= S[end]
    T₀[:, :, z .< depth_of_steps[end]] .= T[end]

    set!(model, S = S₀, T = T₀)

    return nothing
end
function set_staircase_initial_conditions!(model, ics::SmoothStepInitialConditions)

    # TODO: write methods to set smooth changes using the function provied in
    # `ics::SmoothStepInitialConditions`. I can use something like the above
    # normalise then set with the smoothing function
    #set!(model, S = initial_S_steps, T = initial_T_steps)

    return nothing
end
function set_staircase_initial_conditions!(sdns::StaircaseDNS)

    set_staircase_initial_conditions!(sdns.model, sdns.initial_conditions)

    return nothing
end
