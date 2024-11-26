"""
    function set_initial_conditions!(model, ics)
Set initial salinity and temperature staircases in `model`.
"""
function set_initial_conditions!(model, ics::STStaircaseInitialConditions)

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
function set_initial_conditions!(model, ics::SmoothSTStaircaseInitialConditions)

    # TODO: write methods to set smooth changes using the function provied in
    # `ics::SmoothSTStaircaseInitialConditions`. I can use something like the above
    # normalise then set with the smoothing function
    #set!(model, S = initial_S_steps, T = initial_T_steps)

    return nothing
end
set_initial_conditions!(model, ics::Union{SingleInterfaceICs, PeriodoicSingleInterfaceICs}) =
    set_initial_conditions!(model, ics, ics.interface_smoothing)
function set_initial_conditions!(model, ics::SingleInterfaceICs, interface_smoothing::Nothing)

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
set_initial_conditions!(model, ics::PeriodoicSingleInterfaceICs, interface_smoothing::Nothing) = nothing
function set_initial_conditions!(model, ics, interface_smoothing::Type{<:Tanh})

    depth_of_interface = ics.depth_of_interface
    Lz = model.grid.Lz

    S = Array(ics.salinity_values)
    Sᵤ, Sₗ = S
    ΔS = diff(S)[1]
    S₀(x, y, z) = Tanh(Sₗ, ΔS, 100.0, depth_of_interface, abs(Lz))(x, y, z)

    T = Array(ics.temperature_values)
    Tᵤ, Tₗ = T
    ΔT = diff(T)[1]
    T₀(x, y, z) = Tanh(Tₗ, ΔT, 100.0, depth_of_interface, abs(Lz))(x, y, z)

    set!(model, S = S₀, T = T₀)

    return nothing
end

"Fallback --- don't set any noise."
set_noise!(model, noise::Nothing) = nothing

function set_noise!(model, noise::NamedTuple)

    set_noise!(model, noise.tracers)
    set_noise!(model, noise.velocities)

    return nothing
end
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

    return nothing
end
"""
    function set_noise!(model, noise::TracerNoise)
Add initial noise the `tracer` fields. The tracer noise is added to the values already in
the tracer `Field`s so if nothing else is set it will just be noise.
"""
function set_noise!(model, noise::TracerNoise)

    S, T = model.tracers

    S_noise = noise.S_magnitude * randn(size(S))
    S_noise_field = similar(S)
    set!(S_noise_field, S_noise)

    T_noise = noise.T_magnitude * randn(size(T))
    T_noise_field = similar(T)
    set!(T_noise_field, T_noise)

    S₀ = S + S_noise_field
    T₀ = T + T_noise_field

    set!(model, S = S₀, T = T₀)

    return nothing
end
"""
    function set_initial_conditions!(sdns::StaircaseDNS)
Set initial staircase and noise for a `StaircaseDNS`.
"""
function set_initial_conditions!(sdns::StaircaseDNS)

    set_initial_conditions!(sdns.model, sdns.initial_conditions)
    set_noise!(sdns.model, sdns.initial_noise)

    return nothing
end
