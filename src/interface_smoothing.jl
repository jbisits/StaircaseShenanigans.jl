"""
    struct TanhInterfaceSmoothing
Container to set a hyperbolic tangent over a single interface as an initial condition.
"""
mutable struct TanhInterfaceSmoothing{T} <: AbstractInterfaceSmoothing
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
"This is a bit of a hack to allow the steepness to be set without needin to set everything
else which requires knowledge of grid and layer setup. Basically a container for just the
steepness of the tanh change."
Tanh(D) = Tanh(0.0, 0.0, D, 0.0, 0.0)
@inline (p::TanhInterfaceSmoothing)(x, y, z) = p.Cₗ - 0.5 * p.ΔC * (1  + tanh(p.D * (z - p.z_interface) / p.z_range))

Base.summary(p::Type{<:TanhInterfaceSmoothing}) = "tanh smoothing"
Base.summary(p::TanhInterfaceSmoothing) = "tanh smoothing"

struct NoInterfaceSmoothing <: AbstractInterfaceSmoothing end
const NoSmoothing = NoInterfaceSmoothing
Base.summary(ns::Type{<:NoInterfaceSmoothing}) = "no interface smoothing"
