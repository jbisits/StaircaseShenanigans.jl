"Temperature at the _top_ of the domain + ΔT across interface for use as `ValueBoundaryCondition`."
@inline T_jump(i, j, grid, clock, model_fields, ΔT) = @inbounds model_fields.T[i, j, grid.Nz] + ΔT

"Salinity at the _top_ of the domain + ΔS across interface for use as `ValueBoundaryCondition`."
@inline S_jump(i, j, grid, clock, model_fields, ΔS) = @inbounds model_fields.S[i, j, grid.Nz] + ΔS

"Velcotiy at the _top_ of the domain for use as `OpenBoundaryCondition`."
@inline w_jump(i, j, grid, clock, model_fields) = @inbounds model_fields.w[i, j, grid.Nz+1]

"""
    function jump_periodic_boundary_conditions
Add boundary conditions for `S`, `T` and `w` to create *jump periodic* boundary conditions.
The jump peridioc conditions for the tracers is
```math
C(z = Lz) = C(z = 0) + ΔC
```
where ``ΔC`` is the change in tracer across an interface (or multiple interfaces).
For the vertical velocity the boundary condition is
```math
w(z = Lz) = w(z = 0).
```
These boundary conditions aim to provide simulations that do not run down and can be used as
an alternative to triply periodic conditions.
**Note:** These boundary conditions are applied by default in `StaircaseDNS` if the topology
in the z direcion is `Bounded`.
"""
function jump_periodic_boundary_conditions(ics::SingleInterfaceICs, z_topology::Type{<:Bounded})

    ΔT = abs(diff(Array(ics.temperature_values))[1])
    ΔS = abs(diff(Array(ics.salinity_values))[1])
    T_bottom = ValueBoundaryCondition(T_jump, discrete_form=true, parameters = ΔT)
    S_bottom = ValueBoundaryCondition(S_jump, discrete_form=true, parameters = ΔS)
    T_bcs = FieldBoundaryConditions(bottom = T_bottom)
    S_bcs = FieldBoundaryConditions(bottom = S_bottom)

    w_bottom = OpenBoundaryCondition(w_jump, discrete_form=true)
    w_bcs = FieldBoundaryConditions(bottom = w_bottom)

    return (T=T_bcs, S=S_bcs, w=w_bcs)
end
jump_periodic_boundary_conditions(ics::SingleInterfaceICs, z_topology) = nothing
