"""
    function set_staircase_initial_conditions!(model, ics)
Set initial salinity and temperature staircases in `model`.
"""
function set_staircase_initial_conditions!(model, ics::StepInitialConditions)

    depth_of_steps = ics.depth_of_steps
    S = ics.salinity_values
    T = ics.temperature_values
    initial_S_steps(x, y, z) = set_steps(z, S, depth_of_steps)
    initial_T_steps(x, y, z) = set_steps(z, T, depth_of_steps)

    set!(model, S = initial_S_steps, T = initial_T_steps)

    return nothing
end
function set_staircase_initial_conditions!(model, ics::SmoothStepInitialConditions)

    # TODO: write methods to set smooth changes using the function provied in
    # `ics::SmoothStepInitialConditions`.
    #set!(model, S = initial_S_steps, T = initial_T_steps)

    return nothing
end
function set_staircase_initial_conditions!(sdns::StaircaseDNS)

    set_staircase_initial_conditions!(sdns.model, sdns.initial_conditions)

    return nothing
end
