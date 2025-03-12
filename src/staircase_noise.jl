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
"Convenience for all noise at same magnitude `c`, default behaviour is 1e-4"
VelocityNoise(c::Float64=1e-4) = VelocityNoise(c, c, c)
Base.summary(noise::VelocityNoise) = "Random noise magnitude in u, v and w fields $(noise.u_magnitude), $(noise.v_magnitude), and $(noise.w_magnitude) respectively"
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
Base.summary(noise::NamedTuple) = "Random noise in velocity and tracer fields"

const NoiseOptions = Union{<:VelocityNoise{T}, <:TracerNoise{T}} where T
"""
    NoiseAtDepth{D, AN}
Add noise at `depth` or make `depth` an `Array` to set noise over specific region of domain.
"""
struct NoiseAtDepth{A <: AbstractArray, N <: NoiseOptions} <: AbstractNoise
    "Depth range at which to add random noise"
    depth :: A
    "Noise to add either `TracerNoise` or `VelocityNoise`"
    noise :: N
end
Base.summary(nd::NoiseAtDepth{<:AbstractArray, <:TracerNoise}) = "Random noise magnitude in S and T fields $(nd.noise.S_magnitude) and $(nd.noise.T_magnitude) where z ∈ $(nd.depth)"
Base.summary(nd::NoiseAtDepth{<:AbstractArray, <:VelocityNoise}) = "Random noise magnitude in u, v and w fields $(nd.noise.u_magnitude), $(nd.noise.v_magnitude), and $(nd.noise.w_magnitude) where z ∈ $(nd.depth)"

function Base.show(io::IO, noise::AbstractNoise)
    if noise isa VelocityNoise
        println(io, "VelocityNoise")
        println(io, "┣━ u_noise_magnitude: $(noise.u_magnitude)")
        println(io, "┣━ v_noise_magnitude: $(noise.v_magnitude)")
        print(io,   "┗━ w_noise_magnitude: $(noise.w_magnitude)")
    elseif noise isa TracerNoise
        println(io, "TracerNoise")
        println(io, "┣━ S_noise_magnitude: $(noise.S_magnitude)")
        print(io,   "┗━ T_noise_magnitude: $(noise.T_magnitude)")
    elseif noise isa NoiseAtDepth
        println(io, "NoiseAtDepth")
        println(io, "┣━ noise_depth: $(noise.depth)")
        print(io,   "┗━━ noise_type: $(noise.noise)")
    end
end
