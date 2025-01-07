"""
    struct TanhInterfaceSmoothing
Container to set a hyperbolic tangent over a single interface as an initial condition.
"""
struct TanhInterfaceSmoothing{T} <: AbstractInterfaceSmoothing
    "Tracer in the deeper of the two layers seperated by an interface"
    Cₗ :: T
    "Change in tracer content over an interface"
    ΔC :: T
    "Steepness of the tanh change"
     D :: T
    "Location of the interface"
    z_interface :: T
    "Total z range of the two layers seperated by a given interface"
    z_range :: T
end
const Tanh = TanhInterfaceSmoothing{T} where {T}
@inline (p::TanhInterfaceSmoothing)(x, y, z) = p.Cₗ - 0.5 * p.ΔC * (1  + tanh(p.D * (z - p.z_interface) / p.z_range))

"Container for the steepness of the `tanh` over the interface for both salinity and temperature.
Allows easy setting of both and avoids having to set other parts of the container before the
model has been setup."
struct TanhInterfaceSteepness{T}
    "Salinity interface steepness"
    DS :: T
    "Temperature interface steepness"
    DT :: T
end
"Some defaults that are unstable to diffusive convection."
TanhInterfaceSteepness() = TanhInterfaceSteepness(500.0, 250.0)
TanhInterfaceSteepness(S) = TanhInterfaceSteepness(S, S / 3)

Base.summary(p::Type{<:TanhInterfaceSmoothing}) = "tanh smoothing"
Base.summary(p::TanhInterfaceSmoothing) = "tanh smoothing"
Base.summary(p::TanhInterfaceSteepness) = "tanh smoothing"

struct NoInterfaceSmoothing <: AbstractInterfaceSmoothing end
const NoSmoothing = NoInterfaceSmoothing
Base.summary(ns::Type{<:NoInterfaceSmoothing}) = "no interface smoothing"
