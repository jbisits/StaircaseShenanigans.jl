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

**Note:** currently using `model_setup` does not allow boundary conditions to be added.
Instead use `DNSModel` and build a `StaircaseDNS` from the `DNSModel`.
"""
function StaircaseDNS(model_setup::NamedTuple, initial_conditions::SingleInterfaceICs, initial_noise)

    background_fields = S_and_T_background_fields(initial_conditions, model_setup.domain_extent.Lz)
    model = DNSModel(model_setup...; background_fields)

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
- `domain_extent` as a `NamedTuple` in the format `(Lx = , Ly = , Lz = )`;
- `domain_topology`, `NamedTuple` in form `(x = a, y = b, z = c)` where `a`, `b` and `c` are
Oceananigans.jl topologies see [Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/dev/grids/#grids_tutorial)
for more info
- `resolution` as a `NamedTuple` in the format `(Nx = , Ny = , Nz = )`;
- `diffusivities` as a `NamedTuple` in the format `(ν = , κ = )`, **note** to set different
diffusivities for temperature and salinity `κ` must also be a `NamedTuple` in the format
`κ = (S = , T = )`;
- `eos` a `BoussinesqEquationOfState`, default is [`TEOS10EquationOfState`](https://clima.github.io/SeawaterPolynomials.jl/dev/API/#SeawaterPolynomials.TEOS10.TEOS10EquationOfState)
but any of the [`RoquetEquationOfState`s](https://clima.github.io/SeawaterPolynomials.jl/dev/API/#SeawaterPolynomials.SecondOrderSeawaterPolynomials.RoquetEquationOfState)
may be used.


## Keyword arguments:

- `forcing = nothing`, add a restoring value to part of the staircase. Done by passing the
`forcing` argument to `NonhydrostaticModel`. For how to implement `forcing` see the relevant
part of the
[Oceananigans documentation](https://clima.github.io/OceananigansDocumentation/dev/model_setup/forcing_functions/#forcing)
- `zgrid_stretching` stretch the grid in the `z` dimension at the bottom of domain at
the rate `stretching`, `false` by default;
- `refinement = 1.2` spacing near the surface in the `z` dimension;
- `stretching = 100` rate of stretching at the bottom of grid in the `z` dimension.
"""
function DNSModel(architecture, diffusivities::NamedTuple, domain_extent::NamedTuple,
                  domain_topology::NamedTuple, resolution::NamedTuple,
                  eos::BoussinesqEquationOfState=TEOS10EquationOfState();
                  forcing = nothing,
                  boundary_conditions = nothing,
                  background_fields = nothing,
                  zgrid_stretching = false,
                  refinement = 1.05,
                  stretching = 40)

    Lx, Ly, Lz = domain_extent.Lx, domain_extent.Ly, domain_extent.Lz
    x_top, y_top, z_top = domain_topology.x, domain_topology.y, domain_topology.z
    Nx, Ny, Nz = resolution.Nx, resolution.Ny, resolution.Nz
    zgrid = zgrid_stretching ? grid_stretching(-Lz, Nz, refinement, stretching) : (Lz, 0)

    grid = RectilinearGrid(architecture,
                           topology = (x_top, y_top, z_top),
                           size = (Nx, Ny, Nz),
                           x = (-Lx/2, Lx/2),
                           y = (-Ly/2, Ly/2),
                           z = zgrid)

    buoyancy = SeawaterBuoyancy(equation_of_state = eos)
    tracers = (:S, :T)

    closure = ScalarDiffusivity(; ν = diffusivities.ν, κ = diffusivities.κ)

    timestepper = :RungeKutta3

    advection = Centered(order = 2)

    forcing = isnothing(forcing) ? NamedTuple() : forcing

    boundary_conditions = isnothing(boundary_conditions) ? NamedTuple() : boundary_conditions

    background_fields = isnothing(background_fields) ? NamedTuple() : background_fields

    return NonhydrostaticModel(; grid, buoyancy, tracers, closure, timestepper, advection,
                                 forcing, boundary_conditions, background_fields)

end
"""
    function grid_stretching(Lz, Nz; refinement, stretching)
Stretch the vertical coordinate of the grid. This function was taken from an
[Oceananigans.jl example](https://clima.github.io/OceananigansDocumentation/dev/generated/ocean_wind_mixing_and_convection/).
The rate of stretching at the bottom is controlled by the `stretching` argument and the
spacing near the surface is controlled by `refinement`.
"""
function grid_stretching(Lz::Number, Nz::Number, refinement::Number, stretching::Number)

    # Normalise height ranging from 0 to 1
    h(k) = (k - 1) / Nz
    # Linear near-surface generator
    ζ₀(k) = 1 + (h(k) - 1) / refinement
    # Bottom-intensified stretching function
    Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))
    # Generating function
    z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

    return z_faces

end
