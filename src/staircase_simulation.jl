"""
    function SDNS_simulation_setup
Build a `simulation` from `sdns`. Only required arguments are the `sdns` and `stop_time`.
All other `kwargs` have sensible defaults.
"""
function SDNS_simulation_setup(sdns::StaircaseDNS, stop_time::Number,
                                save_custom_output!::Function=no_custom_output!,
                                save_velocities!::Function=no_velocities!;
                                Δt = nothing,
                                max_Δt = Inf,
                                max_change = 1.1,
                                save_schedule = 60, # seconds
                                save_file = :netcdf,
                                output_path = SIMULATION_PATH,
                                checkpointer_time_interval = nothing,
                                cfl = 0.2,
                                diffusive_cfl = 0.75,
                                overwrite_saved_output = true,
                                with_halos = false)

    Δt = isnothing(Δt) ? initial_timestep(sdns) : Δt
    simulation = Simulation(sdns.model; Δt, stop_time)

    # time step adjustments
    wizard = TimeStepWizard(; cfl, diffusive_cfl, max_change, max_Δt)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

    # progress reporting
    simulation.callbacks[:progress] = Callback(simulation_progress, IterationInterval(100))

    output_dir = output_directory(sdns, stop_time, output_path)
    save_info = (save_schedule, save_file, output_dir, overwrite_saved_output, with_halos)

    # model tracers
    save_tracers!(simulation, sdns, save_info)

    # model velocities
    save_velocities!(simulation, sdns, save_info)

    # Custom saved output
    save_custom_output!(simulation, sdns, save_info)

    # Checkpointer setup
    checkpointer_setup!(simulation, sdns, output_dir, checkpointer_time_interval)

    non_dimensional_numbers!(simulation, sdns)

    save_background_state!(simulation, sdns)

    return simulation

end
"""
    function initial_timestep(sdns)
Compute an initial timestep from the diffusivities set in the model. The timestep is
a sensible choice based on the diffusivities set in the simulation.
"""
initial_timestep(sdns) = 0.2 * cell_diffusion_timescale(sdns.model)
function compute_Ra(ics::SingleInterfaceICs, κₛ, κₜ, ν, eos, int_depth, L; g = 9.81)

    salinity = Array(ics.salinity_values)
    temperature = Array(ics.temperature_values)

    return compute_Ra(salinity, temperature, κₛ, κₜ, ν, eos, int_depth, L, g)
end
function compute_Ra(ics::StaircaseICs, κₛ, κₜ, ν, eos, int_depth, L; g = 9.81)

    salinity = Array(ics.salinity_values)
    temperature = Array(ics.temperature_values)
    RaS = Array{Float64}(undef, length(salinity)-1)
    RaT = similar(RaS)

    for i ∈ eachindex(RaS)
        RaS[i], RaT[i] = compute_Ra(salinity[i:i+1], temperature[i:i+1], κₛ, κₜ, ν, eos, int_depth, L, g)
    end

    return RaS, RaT
end
"Salinity and temperature Rayleigh numbers."
function compute_Ra(salinity, temperature, κₛ, κₜ, ν, eos, int_depth, L, g)

    S₀ᵘ, S₀ˡ = salinity
    T₀ᵘ, T₀ˡ = temperature
    S̄ = 0.5 * (S₀ˡ - S₀ᵘ)
    T̄ = 0.5 * (T₀ˡ - T₀ᵘ)
    ΔS = S₀ˡ - S₀ᵘ
    ΔT = T₀ˡ - T₀ᵘ
    S_u = S_g = S₀ᵘ
    S_l = S_f = S₀ˡ
    T_u = T_f = T₀ᵘ
    T_l = T_g = T₀ˡ

    ρ_u = total_density(T_u, S_u, int_depth, eos)
    ρ_l = total_density(T_l, S_l, int_depth, eos)
    ρ_f = total_density(T_f, S_f, int_depth, eos)
    ρ_g = total_density(T_g, S_g, int_depth, eos)

    ρ_m = total_density(T̄, S̄, int_depth, eos)

    # McDougall (1981)
    α̃ = (0.5 * (ρ_f - ρ_g) - 0.5 * (ρ_l - ρ_u)) / (ρ_m * ΔT)
    β̃ = (0.5 * (ρ_f - ρ_g) + 0.5 * (ρ_l - ρ_u)) / (ρ_m * ΔS)

    RaS = (g * β̃ * ΔS * L^3) / (ν * κₛ)
    RaT = (g * α̃ * ΔT * L^3) / (ν * κₜ)

    return RaS, RaT
end
"""
    non_dimensional_numbers!(simulation::Simulation, sdns::StaircaseDNS)
Save these to the output file.
"""
function non_dimensional_numbers!(simulation::Simulation, sdns::StaircaseDNS)

    model, initial_conditions = sdns.model, sdns.initial_conditions
    eos = model.buoyancy.formulation.equation_of_state
    ν = model.closure.ν
    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    Pr = ν / κₜ
    Sc = ν / κₛ
    Le = κₜ / κₛ
    Lz = model.grid.Lz
    int_depth = interface_depth(initial_conditions)
    L = Lz - int_depth #
    RaS, RaT = compute_Ra(initial_conditions, κₛ, κₜ, ν, eos, int_depth, L)

    nd_nums = Dict("Pr" => Pr, "Sc" => Sc, "Le" => Le, "RaS" => RaS, "RaT" => RaT,
                    "Rᵨ" => initial_conditions.R_ρ)

    NCDataset(simulation.output_writers[:computed_output].filepath, "a") do ds
        ds.attrib["interface_depth"] = int_depth
        ds.attrib["EOS"] = summary(eos.seawater_polynomial)
        ds.attrib["Reference density (kgm⁻³)"] = eos.reference_density
        ds.attrib["ν (m²s⁻¹)"]  = ν
        ds.attrib["κₛ (m²s⁻¹)"] = κₛ
        ds.attrib["κₜ (m²s⁻¹)"] = κₜ
        for key ∈ keys(nd_nums)
            ds.attrib[key] = nd_nums[key]
        end
    end

    return nothing

end
interface_depth(ics::SingleInterfaceICs) = ics.depth_of_interface
interface_depth(ics::StaircaseICs) = ics.depth_of_interfaces[end]
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
                       overwrite_saved_output, with_halos)

    model  = sdns.model
    ics = sdns.initial_conditions
    S = !(ics.background_state isa Type{<:NoBackground}) ? model.tracers.S : Field(model.background_fields.tracers.S + model.tracers.S)
    T = !(ics.background_state isa Type{<:NoBackground}) : Field(model.background_fields.tracers.T + model.tracers.T)

    S_ha = Average(S, dims = (1, 2))
    T_ha = Average(T, dims = (1, 2))

    tracers = Dict("S" => S, "S_ha" => S_ha, "T" => T, "T_ha" => T_ha)

    if !(ics.background_state isa Type{<:NoBackground})
        anomalies = Dict("S′" => model.tracers.S, "T′" => model.tracers.T)
        merge!(tracers, anomalies)
    end


    simulation.output_writers[:tracers] =
        save_file == :netcdf ? NetCDFWriter(model, tracers;
                                            filename = "tracers",
                                            dir = output_dir,
                                            overwrite_existing = overwrite_saved_output,
                                            schedule = TimeInterval(save_schedule),
                                            with_halos) :
                                JLD2Writer(model, tracers;
                                           filename = "tracers",
                                           dir = output_dir,
                                           schedule = TimeInterval(save_schedule),
                                           overwrite_existing = overwrite_saved_output,
                                           with_halos)

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
                          overwrite_saved_output, with_halos)

    model = sdns.model
    u, v, w = model.velocities
    velocities = Dict("u" => u, "v" => v, "w" => w)

    simulation.output_writers[:velocities] =
        save_file == :netcdf ? NetCDFWriter(model, velocities;
                                                filename = "velocities",
                                                dir = output_dir,
                                                overwrite_existing = overwrite_saved_output,
                                                schedule = TimeInterval(save_schedule),
                                                with_halos
                                                ) :
                                JLD2Writer(model, velocities;
                                                filename = "velocities",
                                                dir = output_dir,
                                                schedule = TimeInterval(save_schedule),
                                                overwrite_existing = overwrite_saved_output,
                                                with_halos)

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
                                    overwrite_saved_output, with_halos)

    model = sdns.model
    w = model.velocities.w
    velocities = Dict("w" => w)

    simulation.output_writers[:velocities] =
    save_file == :netcdf ? NetCDFWriter(model, velocities;
                                filename = "velocities",
                                dir = output_dir,
                                overwrite_existing = overwrite_saved_output,
                                schedule = TimeInterval(save_schedule),
                                with_halos
                                ) :
                JLD2Writer(model, velocities;
                                filename = "velocities",
                                dir = output_dir,
                                schedule = TimeInterval(save_schedule),
                                overwrite_existing = overwrite_saved_output,
                                with_halos)

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
                               overwrite_saved_output, with_halos, reference_gp_height)

    model = sdns.model
    ics = sdns.initial_conditions
    S = !(ics.background_state isa Type{<:NoBackground}) ? model.tracers.S : Field(model.background_fields.tracers.S + model.tracers.S)
    T = !(ics.background_state isa Type{<:NoBackground}) ? model.tracers.T : Field(model.background_fields.tracers.T + model.tracers.T)

    σ = seawater_density(model, temperature = T, salinity = S, geopotential_height = reference_gp_height)
    σ_ha = Average(σ, dims = (1, 2))
    N² = buoyancy_frequency(model)
    N²_ha = Average(N², dims = (1, 2))
    ε = KineticEnergyDissipationRate(model)
    ∫ε = Integral(ε)
    ε_maximum = Reduction(maximum!, ε, dims = (1, 2, 3))
    Eₖ = KineticEnergy(model)
    ∫Eₖ = Integral(Eₖ)
    wb = BuoyancyProductionTerm(model)
    ∫wb = Integral(wb)

    computed_outputs = Dict("σ" => σ, "N²" => N², "σ_ha" => σ_ha, "N²_ha" => N²_ha,
                            "∫ε" => ∫ε, "ε_maximum" => ε_maximum, "∫Eₖ" => ∫Eₖ, "∫wb" => ∫wb)
    oa = Dict(
        "σ" => Dict("longname" => "Seawater potential density calculated using equation of state in model.",
                    "units" => "kgm⁻³"),
        "N²" => Dict("longname" => "Buoyancy frequency.",
                    "units" => "s⁻¹"),
        "σ_ha" => Dict("longname" => "Horizontal average seawater potential density calculated using equation of state in model.",
                    "units" => "kgm⁻³"),
        "N²_ha" => Dict("longname" => "Horizontal average buoyancy frequency.",
                    "units" => "s⁻¹"),
        "∫ε" => Dict("longname" => "Volume integrated turbulent kinetic energy dissipation rate.",
                        "units" => "m⁵s⁻³"),
        "ε_maximum" => Dict("longname" => "Maximum (in space) local TKE rate.",
                        "units" => "m⁵s⁻³"),
        "∫Eₖ" => Dict("longname" => "Volume integrated turbulent kinetic energy.",
                        "units" => "m⁵s⁻²"),
        "∫wb" => Dict("longname" => "Volume integrated buoyancy flux.",
                        "units" => "m⁵s⁻³")
        )

    if !(ics.background_state isa Type{<:NoBackground})
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
        save_file == :netcdf ? NetCDFWriter(model, computed_outputs;
                                            filename = "computed_output",
                                            dir = output_dir,
                                            overwrite_existing = overwrite_saved_output,
                                            schedule = TimeInterval(save_schedule),
                                            output_attributes = oa,
                                            with_halos) :
                                JLD2Writer(model, computed_outputs;
                                          filename = "computed_output",
                                          dir = output_dir,
                                          schedule = TimeInterval(save_schedule),
                                          overwrite_existing = overwrite_saved_output,
                                          with_halos)


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

    if ow isa NetCDFWriter
        NCDataset(ow.filepath) do ds
            key_present = haskey(ds, key)
        end
    elseif ow isa JLD2Writer
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
