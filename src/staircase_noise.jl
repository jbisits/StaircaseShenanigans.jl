Base.iterate(noise::AbstractNoise, state = 1) =
    state > length(fieldnames(noise)) ? nothing : (getfield(noise, state), state + 1)

"""
    struct VelocityNoise{T}
Container for the magnitudes of normally distributed noise added to velocity fields.
"""
struct VelocityNoise{T} <: AbstractNoise
    "Noise magnitude for `u` velocity field"
    u_magnitude :: T
    "Noise magnitude for `v` velocity field"
    v_magnitude :: T
    "Noise magnitude for `w` velocity field"
    w_magnitude :: T
end
"""
    struct TracerNoise{T}
Container for the magnitudes of normally distributed noise added to tracer fields.
"""
struct TracerNoise{T} <: AbstractNoise
    "Noise magnitude for `S` tracer field"
    S_magnitude :: T
    "Noise magnitude for `T` tracer field"
    T_magnitude :: T
end
"Convenience for all noise at same magnitude `c`, default behaviour is 1e-4"
TracerNoise(c::Float64=1e-4) = TracerNoise(c, c)
Base.summary(noise::TracerNoise) = "Random noise magnitude in S and T fields $(noise.S_magnitude) and $(noise.T_magnitude) respectively"

function Base.show(io::IO, noise::AbstractNoise)
    if noise isa VelocityNoise
        println(io, "VelocityNoise")
        println(io, "┣━━ u_noise_magnitude: $(noise.u_magnitude)")
        println(io, "┣━━ v_noise_magnitude: $(noise.v_magnitude)")
        print(io,   "┗━━ w_noise_magnitude: $(noise.w_magnitude)")
    elseif noise isa TracerNoise
        println(io, "TracerNoise")
        println(io, "┣━━ S_noise_magnitude: $(noise.S_magnitude)")
        print(io,   "┗━━ T_noise_magnitude: $(noise.T_magnitude)")
    end
end
"Convenience for all noise at same magnitude `c`, default behaviour is 1e-4"
VelocityNoise(c::Float64=1e-4) = VelocityNoise(c, c, c)
Base.summary(noise::VelocityNoise) = "Random noise magnitude in u, v and w fields $(noise.u_magnitude), $(noise.v_magnitude), and $(noise.w_magnitude) respectively"
