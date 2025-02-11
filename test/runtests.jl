using StaircaseShenanigans, GibbsSeaWater
using SeawaterPolynomials: thermal_expansion, haline_contraction
using Test

@testset "StaircaseShenanigans.jl" begin

    include("basic_configuation.jl")

    @testset "Setting initial conditions" begin

        include("model_ic_instantiation.jl")
        for interface_smoothing ∈ smoothing
            for background_state ∈ background
                @test SingleInterfaceICs(eos, depth_of_interface, salinity, temperature;
                                        interface_smoothing, background_state) isa STSingleInterfaceInitialConditions
            end
        end

        include("model_ic_setting.jl")
        @test all(unique(interior(sdns_staircase.model.tracers.S, 1, 1, :)) .== reverse(salinity_staircase))
        @test all(unique(interior(sdns_staircase.model.tracers.T, 1, 1, :)) .== reverse(temperature_staircase))
        @test all(unique(interior(sdns_interface.model.tracers.S, 1, 1, :)) .== reverse(salinity))
        @test all(unique(interior(sdns_interface.model.tracers.T, 1, 1, :)) .== reverse(temperature))
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

    @testset "Function for diffusivity closure" begin

        diffusivities = (ν = 1e-5, κ = (S = enhance_κₛ, T = enhance_κₜ),
                        parameters = (diff_change = 0.5, enhance = 10), discrete_form = true)
        model = DNSModel(architecture, diffusivities, domain_extent, domain_topology, resolution)
        stop_time = 5
        simulation = Simulation(model; Δt = 0.1, stop_time)
        @test simulation.model.clock.time > stop_time
    end

end
