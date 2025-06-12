using StaircaseShenanigans, GibbsSeaWater, JLD2
using SeawaterPolynomials: thermal_expansion, haline_contraction
using StaircaseShenanigans: φ_interface_flux!
using Test, CairoMakie

@testset "StaircaseShenanigans.jl" begin

    include("basic_configuation.jl")

    @testset "Setting initial conditions" begin

        include("model_ic_instantiation.jl")
        for interface_smoothing ∈ smoothing
            for background_state ∈ background
                @test SingleInterfaceICs(eos, depth_of_interface, salinity, temperature;
                                         interface_smoothing, background_state) isa STSingleInterfaceInitialConditions
                @test StaircaseICs(eos, number_of_interfaces, depth_of_interfaces, salinities, temperatures;
                                   interface_smoothing, background_state) isa STStaircaseInitialConditions
            end
        end

        include("model_ic_setting.jl")
        @test all(unique(interior(sdns_staircase.model.tracers.S, 1, 1, :)) .== reverse(salinity_staircase))
        @test all(unique(interior(sdns_staircase.model.tracers.T, 1, 1, :)) .== reverse(temperature_staircase))
        @test all(unique(interior(sdns_interface.model.tracers.S, 1, 1, :)) .== reverse(salinity))
        @test all(unique(interior(sdns_interface.model.tracers.T, 1, 1, :)) .== reverse(temperature))

        include("noise_range.jl")
        @test all(interior(sdns_noise_range.model.tracers.S, :, :, z_noise_range) .!= 34.7)
        @test all(interior(sdns_noise_range.model.tracers.S, :, :, z_no_noise_range) .== 34.7)
        @test all(interior(sdns_noise_range.model.tracers.T, :, :, z_noise_range) .!= 0.5)
        @test all(interior(sdns_noise_range.model.tracers.T, :, :, z_no_noise_range) .== 0.5)

        include("background_field_setting.jl")
        S_background_profile = Field(sdns_interface_linear_bf.model.background_fields.tracers.S)
        compute!(S_background_profile)
        @test S_background_profile isa Field
        T_background_profile = Field(sdns_interface_linear_bf.model.background_fields.tracers.T)
        compute!(T_background_profile)
        @test T_background_profile isa Field
        S_background_profile = Field(sdns_interface_tanh_bf.model.background_fields.tracers.S)
        compute!(S_background_profile)
        @test S_background_profile isa Field
        T_background_profile = Field(sdns_interface_tanh_bf.model.background_fields.tracers.T)
        compute!(T_background_profile)
        @test T_background_profile isa Field
    end

    @testset "R_ρ calculation" begin

        include("R_sub_rho.jl")
        ΔS, ΔΘ = diff(salinity[1:2]), diff(temperature[1:2])
        Sₘ, Θₘ = (sum(salinity[1:2]) / 2), (sum(temperature[1:2]) / 2)
        β = haline_contraction(Θₘ, Sₘ, 0, model.buoyancy.formulation.equation_of_state)
        α = -thermal_expansion(Θₘ, Sₘ, 0, model.buoyancy.formulation.equation_of_state)

        @test all(step_ics.R_ρ .≈ (β * ΔS) / (α * ΔΘ))
        @test (interface_ics.R_ρ ≈ vec((β * ΔS) / (α * ΔΘ))[1])

    end

    @testset "CustomLienarEquationOfState initialisation" begin

        for _ in 1:5
            Θ = rand(range(-1.5, 20.5, length = 20))
            S = rand(range(34.5, 34.7, length = 20))
            ρ₀ = 1024.6
            true_coefficients = -ρ₀ * gsw_alpha(S, Θ, 0), ρ₀ * gsw_beta(S, Θ, 0)
            custom_linear = CustomLinearEquationOfState(Θ, S, reference_density = ρ₀)
            sp = custom_linear.seawater_polynomial
            linear_coefficients = sp.R₀₁₀, sp.R₁₀₀
            other_coefficients = sp.R₁₀₁, sp.R₂₀₀, sp.R₀₂₀, sp.R₁₁₀
            @test sp isa SecondOrderSeawaterPolynomial
            @test all(other_coefficients .== 0)
            @test all(linear_coefficients .== true_coefficients)
        end

    end

    @testset "Diagnostics saving and plotting" begin

        ## This does not test accuracy just that a file is saved
        model = DNSModel(architecture, diffusivities, domain_extent, domain_topology, resolution, eos)
        depth_of_interface = -0.5
        salinity = [34.54, 34.70]
        temperature = [-1.5, 0.5]
        interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature,
                                           interface_smoothing = TanhInterfaceThickness())
        sdns = StaircaseDNS(model, interface_ics, nothing)
        set_initial_conditions!(sdns)
        stop_time = 2 * 60 # seconds
        save_schedule = 4  # seconds
        output_path = joinpath(@__DIR__, "output")
        simulation = SDNS_simulation_setup(sdns, stop_time, save_computed_output!, save_vertical_velocities!;
                                            save_schedule,
                                            output_path, max_Δt = 5)
        min_spacing = abs(domain_extent.Lz) / resolution.Nz
        @test simulation.Δt == 0.2 * (min_spacing^2 / sdns.model.closure.ν) # tests `initial_timestep`
        run!(simulation)
        compute_R_ρ!(simulation.output_writers[:computed_output].filepath,
                     simulation.output_writers[:tracers].filepath,
                     (-0.4, -0.2), (-0.8, -0.6), eos)

        # Animations
        animate_density(simulation.output_writers[:computed_output].filepath, "σ", xslice = 2, yslice = 2)
        @test isfile("density_Nsquared.mp4")
        rm("density_Nsquared.mp4")
        animate_tracers(simulation.output_writers[:tracers].filepath, xslice = 2, yslice = 2)
        @test isfile("tracers.mp4")
        rm("tracers.mp4")
        animate_vertical_velocity(simulation.output_writers[:velocities].filepath, xslice = 2, yslice = 2)
        @test isfile("w.mp4")
        rm("w.mp4")

        # Diagnositcs
        diagnostics_file = joinpath(output_path, "test_diagnostics.jld2")
        save_diagnostics!(diagnostics_file,
                          simulation.output_writers[:tracers].filepath,
                          simulation.output_writers[:computed_output].filepath,
                          simulation.output_writers[:velocities].filepath,
                          group = "lineareos")
        @test isfile(diagnostics_file)
        save_diagnostics!(diagnostics_file,
                          simulation.output_writers[:tracers].filepath,
                          simulation.output_writers[:computed_output].filepath,
                          simulation.output_writers[:velocities].filepath,
                          group = "nonlineareos")
        @test isfile(diagnostics_file)

        output = jldopen(diagnostics_file)
        ΔS, ΔT = output["lineareos/ΔS"], output["lineareos/ΔT"]
        Eb, Eb_lower, Eb_upper = output["nonlineareos/Eb"], output["nonlineareos/Eb_lower"], output["nonlineareos/Eb_upper"]
        Ep, Ep_lower, Ep_upper = output["nonlineareos/Ep"], output["nonlineareos/Ep_lower"], output["nonlineareos/Ep_upper"]
        close(output)
        # salinity and temperature averages are correcy
        @test ΔS[1] ≈ salinity[1] - salinity[2]
        @test ΔT[1] ≈ temperature[1] - temperature[2]
        @test ΔS[2] ≈ salinity[1] - salinity[2]
        @test ΔT[2] ≈ temperature[1] - temperature[2]
        # test that the PE's sum to the total
        @test all(Eb .≈ Eb_lower .+ Eb_upper)
        @test all(Ep .≈ Ep_lower .+ Ep_upper)

        update_diagnostic!(diagnostics_file, "nonlineareos", "S_flux",
                          simulation.output_writers[:tracers].filepath,
                          simulation.output_writers[:computed_output].filepath)
        @test isfile(diagnostics_file)

        update_diagnostic!(diagnostics_file, "lineareos", "hₜ",
                            simulation.output_writers[:tracers].filepath,
                            simulation.output_writers[:computed_output].filepath)
        @test isfile(diagnostics_file)

        update_diagnostic!(diagnostics_file, "nonlineareos", "Ẽ",
                            simulation.output_writers[:tracers].filepath,
                            simulation.output_writers[:computed_output].filepath,
                            interface_offset = 2)
        @test isfile(diagnostics_file)

        save_diagnostic!(diagnostics_file, φ_interface_flux!,
                         (simulation.output_writers[:tracers].filepath, :T, "extra_group/"))

        @test isfile(diagnostics_file)

        output = jldopen(diagnostics_file)
        @test all(sort(keys(output)) .== sort(["lineareos", "nonlineareos", "extra_group", "dims"]))
        close(output)

        rm(diagnostics_file)
        rm(simulation.output_writers[:tracers].filepath)
        rm(simulation.output_writers[:computed_output].filepath)
        rm(simulation.output_writers[:velocities].filepath)
    end


    @testset "Function for diffusivity closure" begin

        diffusivities = (ν = 1e-5, κ = (S = enhance_κₛ, T = enhance_κₜ),
                        parameters = (κₛ = 1e-8, κₜ = 1e-6, start_enhance = 1/60, end_enhance = 2/60,
                                      κ_turb = 1e-3),
                        discrete_form = true)
        model = DNSModel(architecture, diffusivities, domain_extent, domain_topology, resolution)
        stop_time = 5
        simulation = Simulation(model; Δt = 0.1, stop_time)
        run!(simulation)
        @test simulation.model.clock.time ≥ stop_time
    end

    @testset "Diffusivity from non-dim numbers" begin

        ν = 1e-6
        diffs = diffusivities_from_ν(ν, Pr = 10)
        @test diffs.ν == ν
        @test diffs.κ.T == 1e-7
        @test diffs.κ.S == 1e-9
        diffs = diffusivities_from_ν(ν, τ = 0.1, Pr = 10)
        @test diffs.κ.T == 1e-7
        @test diffs.κ.S == 1e-8
        @test diffs.ν == ν
        κₜ = 1e-7
        diffs = diffusivities_from_κₜ(κₜ, Pr = 10)
        @test diffs.ν == 1e-6
        @test diffs.κ.S == 1e-9
        @test diffs.κ.T == κₜ
        diffs = diffusivities_from_κₜ(κₜ, τ = 0.1, Pr = 10)
        @test diffs.ν == 1e-6
        @test diffs.κ.S == 1e-8
        @test diffs.κ.T == κₜ
    end
end
