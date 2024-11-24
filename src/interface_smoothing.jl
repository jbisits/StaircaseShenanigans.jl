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

Base.summary(p::Type{<:TanhInterfaceSmoothing}) = "tanh smoothing"

struct NoSmoothing <: AbstractInterfaceSmoothing end
NoSmoothing() = nothing
Base.summary(ns::NoSmoothing) = "no interface smoothing"
