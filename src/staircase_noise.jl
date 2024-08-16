"""
    struct VelocityNoise{T}
Container for magnitudes for normally distributed noise added to velocity fields.
"""
struct VelocityNoise{T} <: AbstractNoise
    "Magnitude for u velocity field"
    u_magnitude :: T
    "Magnitude for v velocity field"
    v_magnitude :: T
    "Magnitude for w velocity field"
    w_magnitude :: T
end
Base.iterate(noise::AbstractNoise, state = 1) =
    state > length(fieldnames(noise)) ? nothing : (getfield(noise, state), state + 1)

function Base.show(io::IO, noise::AbstractNoise)
    if noise isa VelocityNoise
        println(io, "VelocityNoise")
        println(io, "┣━━ u_noise_magnitude: $(noise.u_magnitude)")
        println(io, "┣━━ v_noise_magnitude: $(noise.v_magnitude)")
        print(io,   "┗━━ w_noise_magnitude: $(noise.w_magnitude)")
    end
end
"Convenience for all noise at same magnitude `c`, default behaviour is 1e-4"
VelocityNoise(c::Float64=1e-4) = VelocityNoise(c, c, c)
