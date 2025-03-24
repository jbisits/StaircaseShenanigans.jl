"""
    struct TanhInterfaceSmoothing
Container to set a hyperbolic tangent over a single interface as an initial condition.
"""
struct TanhInterfaceSmoothing{T} <: AbstractInterfaceSmoothing
    "Tracer in the deeper of the two layers seperated by an interface"
    Cₗ :: T
    "Change in tracer content over an interface"
    ΔC :: T
    "Initial interface thickness"
     h :: T
    "Location of the interface"
    z_interface :: T
    "Total z range of the two layers seperated by a given interface"
    z_range :: T
end
const Tanh = TanhInterfaceSmoothing{T} where {T}
@inline (p::TanhInterfaceSmoothing)(x, y, z) = p.Cₗ - 0.5 * p.ΔC * (1  + tanh((z - p.z_interface) / (p.h * p.z_range)))

"Container for the initial thickness of the `tanh` for both salinity and temperature.
Allows easy setting of both and avoids having to set other parts of the container before the
model has been setup."
struct TanhInterfaceThickness{T}
    "Initial salinity interface thickness"
    hₛ :: T
    "Initial salinity interface thickness"
    hₜ :: T
end
"Some defaults that are unstable to diffusive convection."
TanhInterfaceThickness() = TanhInterfaceThickness(0.05, 3*0.05)
TanhInterfaceThickness(hₛ) = TanhInterfaceThickness(hₛ, 3*hₛ)

Base.summary(p::Type{<:TanhInterfaceSmoothing}) = "tanh smoothing"
Base.summary(p::TanhInterfaceSmoothing) = "tanh smoothing"
Base.summary(p::TanhInterfaceThickness) = "tanh smoothing"

struct NoInterfaceSmoothing <: AbstractInterfaceSmoothing end
const NoSmoothing = NoInterfaceSmoothing
Base.summary(ns::Type{<:NoInterfaceSmoothing}) = "no interface smoothing"
