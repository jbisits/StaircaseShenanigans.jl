"""
    function set_staircase_initial_conditions!(model, ics)
Set initial salinity and temperature staircases in `model`.
"""
function set_staircase_initial_conditions!(model, ics::STStaircaseInitialConditions)

    depth_of_interfaces = ics.depth_of_interfaces
    z = znodes(model.grid, Center())
    S = ics.salinity_values
    T = ics.temperature_values

    S₀ = similar(interior(model.tracers.S, :, :, :))
    T₀ = similar(interior(model.tracers.T, :, :, :))
    S₀[:, :, z .> depth_of_interfaces[1]] .= S[1]
    T₀[:, :, z .> depth_of_interfaces[1]] .= T[1]
    for i ∈ 1:length(depth_of_interfaces)-1
        S₀[:, :, depth_of_interfaces[i+1] .≤ z .< depth_of_interfaces[i]] .= S[i+1]
        T₀[:, :, depth_of_interfaces[i+1] .≤ z .< depth_of_interfaces[i]] .= T[i+1]
    end
    S₀[:, :, z .< depth_of_interfaces[end]] .= S[end]
    T₀[:, :, z .< depth_of_interfaces[end]] .= T[end]

    set!(model, S = S₀, T = T₀)

    return nothing
end
function set_staircase_initial_conditions!(model, ics::SingleInterfaceICs)

    depth_of_interface = ics.depth_of_interface
    z = znodes(model.grid, Center())
    S = ics.salinity_values
    T = ics.temperature_values

    S₀ = similar(interior(model.tracers.S, :, :, :))
    T₀ = similar(interior(model.tracers.T, :, :, :))
    S₀[:, :, z .> depth_of_interface[1]] .= S[1]
    T₀[:, :, z .> depth_of_interface[1]] .= T[1]
    S₀[:, :, z .< depth_of_interface[end]] .= S[2]
    T₀[:, :, z .< depth_of_interface[end]] .= T[2]

    set!(model, S = S₀, T = T₀)

    return nothing
end
"Here the `BackgroundField` behaves as the `initial_condition` and noise is added to the
tracer fields to create an instability."
set_staircase_initial_conditions!(model, ics::PeriodoicSingleInterfaceICs) = nothing
function set_staircase_initial_conditions!(model, ics::SmoothSTStaircaseInitialConditions)

    # TODO: write methods to set smooth changes using the function provied in
    # `ics::SmoothSTStaircaseInitialConditions`. I can use something like the above
    # normalise then set with the smoothing function
    #set!(model, S = initial_S_steps, T = initial_T_steps)

    return nothing
end
set_noise!(model, noise::Nothing) = nothing
"""
    function set_noise!(model, noise::VelocityNoise)
Add initial noise the `velocity` fields.
"""
function set_noise!(model, noise::VelocityNoise)

    u, v, w = model.velocities
    u_noise = noise.u_magnitude * randn(size(u))
    v_noise = noise.v_magnitude * randn(size(v))
    w_noise = noise.w_magnitude * randn(size(w))

    set!(model, u = u_noise, v = v_noise, w = w_noise)

end
"""
    function set_noise!(model, noise::TracerNoise)
Add initial noise the `tracer` fields.
"""
function set_noise!(model, noise::TracerNoise)

    S, T = model.tracers
    S_noise = noise.S_magnitude * randn(size(S))
    T_noise = noise.T_magnitude * randn(size(T))

    set!(model, S = S_noise, T = T_noise)

end
"""
    function set_staircase_initial_conditions!(sdns::StaircaseDNS)
Set initial staircase and noise for a `StaircaseDNS`.
"""
function set_staircase_initial_conditions!(sdns::StaircaseDNS)

    set_staircase_initial_conditions!(sdns.model, sdns.initial_conditions)
    set_noise!(sdns.model, sdns.initial_noise)

    return nothing
end
