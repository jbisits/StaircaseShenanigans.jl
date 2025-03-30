struct StaircaseDNS{NHM <: NonhydrostaticModel,
                    SIC <: AbstractInitialConditions,
                     AN <: Union{AbstractNoise, Nothing, NamedTuple}} <: AbstractStaircaseModel
    "An [Oceananigans.jl `NonhydrostaticModel`](https://clima.github.io/OceananigansDocumentation/dev/appendix/library/#Oceananigans.Models.NonhydrostaticModels.NonhydrostaticModel-Tuple{})"
    model :: NHM
    "Initial conditions"
    initial_conditions :: SIC
    "Initial noise to create instability"
    initial_noise :: AN
end
function Base.show(io::IO, sdns::StaircaseDNS)
    println(io, "StaircaseDirectNumericalSimulation")
    println(io, "┣━━━━━━━━━━━━━━ model: $(summary(sdns.model))")
    println(io, "┣━ initial_conditions: $(summary(sdns.initial_conditions))")
    print(io,   "┗━━━━━━ initial_noise: $(summary(sdns.initial_noise))")
end
"""
    function StaircaseDNS(model, initial_conditions; initial_noise = nothing)
Initialise by passing a `model` that has already been built and set the initial conditions.
"""
function StaircaseDNS(model, initial_conditions; initial_noise = nothing)

    sdns = StaircaseDNS(model, initial_conditions, initial_noise)
    set_initial_conditions!(sdns)

    return sdns
end

"""
    function StaircaseDNS(model_setup::NamedTuple, initial_conditions, initial_noise)
Build the model from `model_setup` then return a `StaircaseDNS`. the `NamedTuple` `model_setup`
requires:
- `architecture` either `CPU` or `GPU`
- `diffusivities`, a `NamedTuple` in the form `(ν = val, κ = (S = val, T = val))` where `val`
is the diffusivity values
- `domain_extent`, a `NamedTuple` in form `(Lx = a, Ly = b, Lz = c)` where `a`, `b` and `c` are
the lengths of the domain
- `domain_topology`, `NamedTuple` in form `(x = a, y = b, z = c)` where `a`, `b` and `c` are
Oceananigans.jl topologies see [Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/dev/grids/#grids_tutorial)
for more info
- `resolution`, a `NamedTuple` in form `(Nx = a, Ny = b, Nz = c)` where `a`, `b` and `c` are
the values of resolution along each dimension
- `eos`, a `BoussinesqEquationOfState` from [SeawaterPolynomials.jl](https://github.com/CliMA/SeawaterPolynomials.jl)

**Note:** By default if the domain is `Bounded` in the z-direction `jump_periodic_boundary_conditions`
are applied.
"""
function StaircaseDNS(model_setup::NamedTuple, initial_conditions::SingleInterfaceICs, initial_noise)

    Lz = model_setup.domain_extent.Lz
    z_range = Lz isa Tuple ? Lz : (Lz, 0)
    background_fields = S_and_T_background_fields(initial_conditions, z_range)
    boundary_conditions = jump_periodic_boundary_conditions(initial_conditions, model_setup.domain_topology.z)
    model = DNSModel(model_setup...; background_fields, boundary_conditions)

    sdns = StaircaseDNS(model, initial_conditions, initial_noise)
    set_initial_conditions!(sdns)

    return sdns
end
# TODO: create methods that use either triple periodic or some kind of restoring for
# `STStaircaseInitialConditions`.
Base.iterate(sdns::AbstractStaircaseModel, state = 1) =
    state > length(fieldnames(sdns)) ? nothing : (getfield(sdns, state), state + 1)

"""
    function DNSModel(architecture, domain_extent::NamedTuple, resolution::NamedTuple,
                      diffusivities::NamedTuple)
Setup a Direct Numerical Simulation [`Model`](https://clima.github.io/OceananigansDocumentation/dev/appendix/library/#Oceananigans.Models.NonhydrostaticModels.NonhydrostaticModel-Tuple{})
on `architecture` (`CPU()` or `GPU()`) over the
`domain_extent` with `resolution` and scalar diffusivities for momentum
(kinematic viscosity `ν`) and the temperature and salinity tracers (`κₜ` and `κₛ`).
The salinity and temperature tracerse will then be evolved with the equation of state `eos`
chosen.

## Function arguments:

- `architecture`, `CPU()` or `GPU()`;
- `domain_extent` as a `NamedTuple` in the format `(Lx = , Ly = , Lz = )`, x and y directions
are setup as periodic and z direction can be either a `Tuple` with ``z ∈ [Lz[1], Lz[2]]`` or
``z ∈ [Lz, 0]`` depending on type passed to `domain_extent`.;
- `domain_topology`, `NamedTuple` in form `(x = a, y = b, z = c)` where `a`, `b` and `c` are
Oceananigans.jl topologies see [Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/dev/grids/#grids_tutorial)
for more info
- `resolution` as a `NamedTuple` in the format `(Nx = , Ny = , Nz = )`;
- `diffusivities` as a `NamedTuple` in the format `(ν = , κ = )`, **note** to set different
diffusivities for temperature and salinity `κ` must also be a `NamedTuple` in the format
`κ = (S = , T = )`. Functions can also be passed but `discrete_form` and `parameters` need to
be added to the `NamedTuple`;
- `eos` a `BoussinesqEquationOfState`, default is [`TEOS10EquationOfState`](https://clima.github.io/SeawaterPolynomials.jl/dev/API/#SeawaterPolynomials.TEOS10.TEOS10EquationOfState)
but any of the [`RoquetEquationOfState`s](https://clima.github.io/SeawaterPolynomials.jl/dev/API/#SeawaterPolynomials.SecondOrderSeawaterPolynomials.RoquetEquationOfState)
may be used.


## Keyword arguments:

- `forcing = nothing`, add a restoring value to part of the staircase. Done by passing the
`forcing` argument to `NonhydrostaticModel`. For how to implement `forcing` see the relevant
part of the
[Oceananigans documentation](https://clima.github.io/OceananigansDocumentation/dev/model_setup/forcing_functions/#forcing)
"""
function DNSModel(architecture, diffusivities::NamedTuple, domain_extent::NamedTuple,
                  domain_topology::NamedTuple, resolution::NamedTuple,
                  eos::BoussinesqEquationOfState=TEOS10EquationOfState();
                  forcing = nothing,
                  boundary_conditions = nothing,
                  background_fields = nothing,
                  centred_advection_scheme_order = 2)

    Lx, Ly, Lz = domain_extent.Lx, domain_extent.Ly, domain_extent.Lz
    x_top, y_top, z_top = domain_topology.x, domain_topology.y, domain_topology.z
    loc = z_top isa Periodic ? (Center(), Center(), Center()) : (Center(), Center(), Face())
    Nx, Ny, Nz = resolution.Nx, resolution.Ny, resolution.Nz
    zgrid = Lz isa Tuple ? Lz : (Lz, 0)

    grid = RectilinearGrid(architecture,
                           topology = (x_top, y_top, z_top),
                           size = (Nx, Ny, Nz),
                           x = (-Lx/2, Lx/2),
                           y = (-Ly/2, Ly/2),
                           z = zgrid)

    buoyancy = SeawaterBuoyancy(equation_of_state = eos)
    tracers = (:S, :T)

    closure = length(diffusivities) > 2 ? ScalarDiffusivity(; ν = diffusivities.ν, κ = diffusivities.κ,
                                                            discrete_form = diffusivities.discrete_form,
                                                            loc,
                                                            parameters = diffusivities.parameters) :
                                          ScalarDiffusivity(ν = diffusivities.ν, κ = diffusivities.κ)

    timestepper = :RungeKutta3

    advection = Centered(order = centred_advection_scheme_order)

    forcing = isnothing(forcing) ? NamedTuple() : forcing

    boundary_conditions = isnothing(boundary_conditions) ? NamedTuple() : boundary_conditions

    background_fields = isnothing(background_fields) ? NamedTuple() : background_fields

    return NonhydrostaticModel(; grid, buoyancy, tracers, closure, timestepper, advection,
                                 forcing, boundary_conditions, background_fields)

end
"""
    enhance_κₛ(i, j, k, grid, clock, fields, p)
Enhance the isotropic, scalar diffusiviy applied to the salintiy field by `p.enhance * 100`,
to give same enhancement to enhanced diffusivity as temperature, after time = `p.diff_change`.
Note `p.diff_change` is expeceted in minutes.
"""
@inline enhance_κₛ(i, j, k, grid, clock, fields, p) = clock.time < p.diff_change * 60 ? p.κₛ : p.κₛ * p.enhance * 100
"""
    enhance_enhance_κₜ(i, j, k, grid, clock, fields, p)
Enhance the isotropic, scalar diffusiviy applied to the temperature field by `p.enhance` after
time `p.diff_change`. Note `p.diff_change` is expeceted in minutes.
"""
@inline enhance_κₜ(i, j, k, grid, clock, fields, p) = clock.time < p.diff_change * 60 ? p.κₜ : p.κₜ * p.enhance
