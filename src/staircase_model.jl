struct StaircaseDNS{NHM <: NonhydrostaticModel,
                    SIC <: AbstractInitialConditions,
                     AN <: Union{AbstractNoise, Nothing, Tuple{AbstractNoise, AbstractNoise}}} <: AbstractStaircaseModel
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
Build the model from `model_setup` then return a `StaircaseDNS`, mainly used to build
the `model` with [reentrant_boundary_conditions](@ref) `boundary_conditions` from `initial_conditions`.
The initial conditions are set after building the `model`.
"""
function StaircaseDNS(model_setup::NamedTuple, initial_conditions::SingleInterfaceICs, initial_noise)

    architecture, diffusivities, domain_extent, resolution, eos =  model_setup
    model = DNSModel(architecture, domain_extent, resolution, diffusivities, eos)

    sdns = StaircaseDNS(model, initial_conditions, initial_noise)
    set_initial_conditions!(sdns)

    return sdns
end

"""
    function StaircaseDNS(model_setup::NamedTuple, initial_conditions, initial_noise)
`StaircaseDNS` that is `Periodic` in z direction and evoloves an anomaly around a
background state. The initial conditions are set after building the `model`.
"""
function StaircaseDNS(model_setup::NamedTuple, initial_conditions::PeriodoicSingleInterfaceICs, initial_noise)

    architecture, diffusivities, domain_extent, resolution, eos = model_setup
    z_topology = Periodic
    background_fields = isnothing(initial_conditions.background_state) ? NamedTuple() :
                            S_and_T_background_fields(initial_conditions, domain_extent.Lz)

    model = DNSModel(architecture, domain_extent, resolution, diffusivities, eos;
                     z_topology, background_fields)

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
function DNSModel(architecture, domain_extent::NamedTuple, resolution::NamedTuple,
                  diffusivities::NamedTuple, eos::BoussinesqEquationOfState=TEOS10EquationOfState();
                  forcing = nothing,
                  boundary_conditions = nothing,
                  background_fields = nothing,
                  z_topology = Bounded,
                  zgrid_stretching = false,
                  refinement = 1.05,
                  stretching = 40)

    Lx, Ly, Lz = domain_extent.Lx, domain_extent.Ly, domain_extent.Lz
    Nx, Ny, Nz = resolution.Nx, resolution.Ny, resolution.Nz
    zgrid = zgrid_stretching ? grid_stretching(-Lz, Nz, refinement, stretching) : (Lz, 0)

    grid = RectilinearGrid(architecture,
                           topology = (Periodic, Periodic, z_topology),
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
"""
    function SDNS_simulation_setup
Build a `simulation` from `sdns`. Only required arguments are the `sdns` and `stop_time`.
All other `kwargs` have defaults (that can be changed).
"""
function SDNS_simulation_setup(sdns::StaircaseDNS, stop_time::Number,
                                save_custom_output!::Function=no_custom_output!,
                                save_velocities!::Function=no_velocities!;
                                Δt = 1e-3,
                                max_Δt = 1e-1,
                                max_change = 1.2,
                                save_schedule = 60, # seconds
                                save_file = :netcdf,
                                output_path = SIMULATION_PATH,
                                checkpointer_time_interval = nothing,
                                cfl = 0.75,
                                diffusive_cfl = 0.75,
                                overwrite_saved_output = true)

    simulation = Simulation(sdns.model; Δt, stop_time)

    # time step adjustments
    wizard = TimeStepWizard(; cfl, diffusive_cfl, max_change, max_Δt)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

    # progress reporting
    simulation.callbacks[:progress] = Callback(simulation_progress, IterationInterval(100))

    output_dir = output_directory(sdns, stop_time, output_path)
    save_info = (save_schedule, save_file, output_dir, overwrite_saved_output)

    # model tracers
    save_tracers!(simulation, sdns, save_info)

    # model velocities
    save_velocities!(simulation, sdns, save_info)

    # Custom saved output
    save_custom_output!(simulation, sdns, save_info)

    # Checkpointer setup
    checkpointer_setup!(simulation, sdns, output_dir, checkpointer_time_interval)

    save_R_ρ!(simulation, sdns)

    save_background_state!(simulation, sdns)

    return simulation

end
"""
    save_R_ρ!(simulation::Simulation, sdns::StaircaseDNS)
To the output file.
"""
function save_R_ρ!(simulation::Simulation, sdns::StaircaseDNS)

    if simulation.output_writers[:tracers] isa NetCDFOutputWriter
        NCDataset(simulation.output_writers[:computed_output].filepath, "a") do ds
            ds.attrib["Step interface R_ρ"] = sdns.initial_conditions.R_ρ
        end
    elseif simulation.output_writers[:tracers] isa JLD2OutputWriter
        jldopen(simulation.output_writers[:computed_output].filepath, "a+") do f
            f["Step interface R_ρ"] = sdns.initial_conditions.R_ρ
        end
    end

    return nothing

end
"""
    function output_directory(sdns::StaircaseDNS, stop_time::Number, output_path)
Create an `output_directory` for saved output based on the `initial_conditions` and
`tracer_perturbation` and length of the simulation.
"""
function output_directory(sdns::StaircaseDNS, stop_time::Number, output_path)

    ic_type = typeof(sdns.initial_conditions)
    ic_string = ic_type <: STStaircaseInitialConditions ? "staircase" : "single_interface"

    eos_string = is_linear_eos(sdns.model.buoyancy.model.equation_of_state.seawater_polynomial)

    stop_time_min = stop_time / 60 ≥ 1 ? string(round(Int, stop_time / 60)) :
                                         string(round(stop_time / 60; digits = 2))
    expt_dir = eos_string *"_"* ic_string *"_"* stop_time_min * "min"

    output_dir = mkpath(joinpath(output_path, expt_dir))
    # make a simulation directory if one is not present
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    return output_dir

end
is_linear_eos(eos::TEOS10SeawaterPolynomial) = "nonlineareos"
function is_linear_eos(eos::SecondOrderSeawaterPolynomial)

    non_linear_coefficients = (eos.R₀₁₁, eos.R₀₂₀, eos.R₁₀₁, eos.R₁₁₀, eos.R₂₀₀)
    eos_type = all(non_linear_coefficients .== 0) ? "lineareos" : "nonlineareos"

    return eos_type
end
"""
    function save_tracers!(simulation, sdns, save_schedule, save_file, output_dir, overwrite_saved_output)
Save `model.tracers` during a `Simulation` using an `OutputWriter`.
"""
function save_tracers!(simulation, sdns, save_schedule, save_file, output_dir,
                       overwrite_saved_output)

    model  = sdns.model
    ics = sdns.initial_conditions
    S = ics isa PeriodoicSingleInterfaceICs ? Field(model.background_fields.tracers.S + model.tracers.S) : model.tracers.S
    T = ics isa PeriodoicSingleInterfaceICs ? Field(model.background_fields.tracers.T + model.tracers.T) : model.tracers.T

    Sᵤ_mean = Average(condition_operand(identity, S, upper_quarter, 0))
    Sₗ_mean = Average(condition_operand(identity, S, lower_quarter, 0))
    Tᵤ_mean = Average(condition_operand(identity, T, upper_quarter, 0))
    Tₗ_mean = Average(condition_operand(identity, T, lower_quarter, 0))

    tracers = Dict("S" => S, "Sᵤ_mean" => Sᵤ_mean, "Sₗ_mean" => Sₗ_mean,
                   "T" => T, "Tᵤ_mean" => Tᵤ_mean, "Tₗ_mean" => Tₗ_mean)

    if ics isa PeriodoicSingleInterfaceICs
        anomalies = Dict("S′" => model.tracers.S, "T′" => model.tracers.T)
        merge!(tracers, anomalies)
    end


    simulation.output_writers[:tracers] =
        save_file == :netcdf ? NetCDFOutputWriter(model, tracers;
                                                  filename = "tracers",
                                                  dir = output_dir,
                                                  overwrite_existing = overwrite_saved_output,
                                                  schedule = TimeInterval(save_schedule)
                                                  ) :
                                JLD2OutputWriter(model, tracers;
                                                 filename = "tracers",
                                                 dir = output_dir,
                                                 schedule = TimeInterval(save_schedule),
                                                 overwrite_existing = overwrite_saved_output)

    return nothing

end
save_tracers!(simulation, sdns, save_info::Tuple) =
    save_tracers!(simulation, sdns, save_info...)
"""
    function save_all_velocities!(simulation, model, save_schedule, save_file, output_dir,
                              overwrite_saved_output)
Save `model.velocities` during a `Simulation` using an `OutputWriter`.
"""
function save_all_velocities!(simulation, sdns, save_schedule, save_file, output_dir,
                          overwrite_saved_output)

    model = sdns.model
    u, v, w = model.velocities
    velocities = Dict("u" => u, "v" => v, "w" => w)

    simulation.output_writers[:velocities] =
        save_file == :netcdf ? NetCDFOutputWriter(model, velocities;
                                                filename = "velocities",
                                                dir = output_dir,
                                                overwrite_existing = overwrite_saved_output,
                                                schedule = TimeInterval(save_schedule)
                                                ) :
                                JLD2OutputWriter(model, velocities;
                                                filename = "velocities",
                                                dir = output_dir,
                                                schedule = TimeInterval(save_schedule),
                                                overwrite_existing = overwrite_saved_output)

    return nothing

end
save_all_velocities!(simulation, sdns, save_info::Tuple) =
    save_all_velocities!(simulation, sdns, save_info...)
"""
    function save_vertical_velocities!(simulation, sdns, save_schedule, save_file, output_dir,
                                    overwrite_saved_output)
Only save vertical velocity.
"""
function save_vertical_velocities!(simulation, sdns, save_schedule, save_file, output_dir,
                                    overwrite_saved_output)

    model = sdns.model
    w = model.velocities.w
    velocities = Dict("w" => w)

    simulation.output_writers[:velocities] =
    save_file == :netcdf ? NetCDFOutputWriter(model, velocities;
                                filename = "velocities",
                                dir = output_dir,
                                overwrite_existing = overwrite_saved_output,
                                schedule = TimeInterval(save_schedule)
                                ) :
                JLD2OutputWriter(model, velocities;
                                filename = "velocities",
                                dir = output_dir,
                                schedule = TimeInterval(save_schedule),
                                overwrite_existing = overwrite_saved_output)

    return nothing

end
save_vertical_velocities!(simulation, sdns, save_info::Tuple) =
    save_vertical_velocities!(simulation, sdns, save_info...)
"Default"
no_velocities!(simulation, sdns, save_info...) = nothing
"""
    function save_computed_output!(simulation, sdns, save_schedule, save_file, output_dir,
                                   overwrite_saved_output, reference_gp_height)
Save selection of computed output during a `Simulation` using an `OutputWriter`.
"""
function save_computed_output!(simulation, sdns, save_schedule, save_file, output_dir,
                               overwrite_saved_output, reference_gp_height)

    model = sdns.model
    ics = sdns.initial_conditions
    S = ics isa PeriodoicSingleInterfaceICs ? Field(model.background_fields.tracers.S + model.tracers.S) : model.tracers.S
    T = ics isa PeriodoicSingleInterfaceICs ? Field(model.background_fields.tracers.T + model.tracers.T) : model.tracers.T

    σ = seawater_density(model, temperature = T, salinity = S, geopotential_height = reference_gp_height)

    computed_outputs = Dict("σ" => σ)
    oa = Dict(
        "σ" => Dict("longname" => "Seawater potential density calculated using equation of state in model.",
                    "units" => "kgm⁻³")
        )

    simulation.output_writers[:computed_output] =
        save_file == :netcdf ? NetCDFOutputWriter(model, computed_outputs;
                                                filename = "computed_output",
                                                dir = output_dir,
                                                overwrite_existing = overwrite_saved_output,
                                                schedule = TimeInterval(save_schedule),
                                                output_attributes = oa
                                                ) :
                                JLD2OutputWriter(model, computed_outputs;
                                                filename = "computed_output",
                                                dir = output_dir,
                                                schedule = TimeInterval(save_schedule),
                                                overwrite_existing = overwrite_saved_output)


    return nothing

end
save_computed_output!(simulation, sdns, save_info::Tuple; reference_gp_height = 0) =
    save_computed_output!(simulation, sdns, save_info..., reference_gp_height)
"Default function for `save_custom_output!` in `sdns_simulation_setup`."
no_custom_output!(simulation, sdns, save_info...) = nothing
"""
    function checkpointer_setup!(simulation, model, output_path, checkpointer_time_interval)
Setup a `Checkpointer` at `checkpointer_time_interval` for a `simulation`
"""
function checkpointer_setup!(simulation, sdns, output_dir,
                             checkpointer_time_interval::Number)

    checkpoint_dir = joinpath(output_dir, "model_checkpoints/")
    isdir(checkpoint_dir) ? nothing : mkdir(checkpoint_dir)
    schedule = TimeInterval(checkpointer_time_interval)
    cleanup = true
    checkpointer = Checkpointer(sdns.model; schedule, dir = checkpoint_dir, cleanup)
    simulation.output_writers[:checkpointer] = checkpointer

    return nothing

end
checkpointer_setup!(simulation, sdns, output_dir, checkpointer_time_interval::Nothing) = nothing
"""
    function simulation_progress(sim)
Useful progress messaging for simulation runs. This function is from an
[Oceananigans.jl example](https://clima.github.io/OceananigansDocumentation/dev/generated/horizontal_convection/#A-progress-messenger).
"""
simulation_progress(sim) = @printf("i: % 6d, sim time: % 1.3f, wall time: % 10s, Δt: % 1.4f, advective CFL: %.2e, diffusive CFL: %.2e\n",
                                    iteration(sim), time(sim), prettytime(sim.run_wall_time),
                                    sim.Δt, AdvectiveCFL(sim.Δt)(sim.model),
                                    DiffusiveCFL(sim.Δt)(sim.model))
