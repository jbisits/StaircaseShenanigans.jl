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
                                cfl = 0.2,
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

    eos_string = is_linear_eos(sdns.model.buoyancy.formulation.equation_of_state.seawater_polynomial)

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
    S = typeof(ics.background_state) <: NoBackground ? model.tracers.S : Field(model.background_fields.tracers.S + model.tracers.S)
    T = typeof(ics.background_state) <: NoBackground ? model.tracers.T : Field(model.background_fields.tracers.T + model.tracers.T)

    Sᵤ_mean = Average(condition_operand(identity, S, upper_quarter, 0))
    Sₗ_mean = Average(condition_operand(identity, S, lower_quarter, 0))
    Tᵤ_mean = Average(condition_operand(identity, T, upper_quarter, 0))
    Tₗ_mean = Average(condition_operand(identity, T, lower_quarter, 0))

    tracers = Dict("S" => S, "Sᵤ_mean" => Sᵤ_mean, "Sₗ_mean" => Sₗ_mean,
                   "T" => T, "Tᵤ_mean" => Tᵤ_mean, "Tₗ_mean" => Tₗ_mean)

    if !(typeof(ics.background_state) <: NoBackground)
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
    S = typeof(ics.background_state) <: NoBackground ? model.tracers.S : Field(model.background_fields.tracers.S + model.tracers.S)
    T = typeof(ics.background_state) <: NoBackground ? model.tracers.T : Field(model.background_fields.tracers.T + model.tracers.T)

    σ = seawater_density(model, temperature = T, salinity = S, geopotential_height = reference_gp_height)
    N² = buoyancy_frequency(model)

    computed_outputs = Dict("σ" => σ, "N²" => N²)
    oa = Dict(
        "σ" => Dict("longname" => "Seawater potential density calculated using equation of state in model.",
                    "units" => "kgm⁻³"),
        "N²" => Dict("longname" => "Buoyancy frequency.",
                    "units" => "s⁻¹")
        )

    if !(typeof(ics.background_state) <: NoBackground)
        S′ = model.tracers.S
        T′ = model.tracers.T
        σ′ = seawater_density(model, temperature = T′, salinity = S′, geopotential_height = reference_gp_height)
        density_anomaly = Dict("σ′" => σ′)
        merge(computed_outputs, density_anomaly)
        oa_ = Dict(
            "σ′" => Dict("longname" => "Seawater potential density calculated from salintiy and temperature anomalies equation of state in model.",
                        "units" => "kgm⁻³")
            )
        merge(oa, oa_)
    end
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
    function check_key_present(simulation, name::Symbol, key)
Check if `key` is already present in `simlulation.output_writers[name]`. Return `true` if
`key` is present and `false` otherwise.
"""
function check_key_present(simulation, name::Symbol, key)

    ow = simulation.output_writers[name]
    key_present = false

    if ow isa NetCDFOutputWriter
        NCDataset(ow.filepath) do ds
            key_present = haskey(ds, key)
        end
    elseif ow isa JLD2OutputWriter
        jldopen(ow.filepath, "a+") do f
            key_present = haskey(f, key)
        end
    end

    return key_present
end
"""
    function simulation_progress(sim)
Useful progress messaging for simulation runs. This function is from an
[Oceananigans.jl example](https://clima.github.io/OceananigansDocumentation/dev/generated/horizontal_convection/#A-progress-messenger).
"""
simulation_progress(sim) = @printf("i: % 6d, sim time: % 1.3f, wall time: % 10s, Δt: % 1.4f, advective CFL: %.2e, diffusive CFL: %.2e\n",
                                    iteration(sim), time(sim), prettytime(sim.run_wall_time),
                                    sim.Δt, AdvectiveCFL(sim.Δt)(sim.model),
                                    DiffusiveCFL(sim.Δt)(sim.model))
