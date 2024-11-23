
"""
    struct TanhInterfaceSmoothing
Container to set a hyperbolic tangent over a single interface as an initial condition.
"""
struct TanhInterfaceSmoothing{T, A} <: AbstractInterfaceSmoothing
    "Tracer in the deeper of the two layers seperated by an interface"
    Cₗ :: T
    "Change in tracer content over an interface"
    ΔC :: T
    "Steepness of the tanh change"
     D :: T
    "Location of the interface"
    z_interface :: T
    "Total z range of the two layers seperated by a given interface"
    z_range :: A
end
@inline (p::TanhInterfaceSmoothing)(x, y, z) = p.Cₗ - 0.5 * p.ΔC * (1  + tanh(p.D * (z - p.z_interface) / p.Lz))

Base.summary(p::TanhInterfaceSmoothing) = "tanh smoothing."
Base.summary(p::Nothing) = "No smoothing."
