using StaircaseShenanigans, GibbsSeaWater
using SeawaterPolynomials: thermal_expansion, haline_contraction
using Test

@testset "StaircaseShenanigans.jl" begin

    include("basic_configuation.jl")

    @testset "Setting initial conditions" begin

        include("model_ic_setup.jl")
        @test all(unique(interior(sdns.model.tracers.S, 1, 1, :)) .== reverse(salinity))
        @test all(unique(interior(sdns.model.tracers.T, 1, 1, :)) .== reverse(temperature))

    end

    @testset "R_ρ calculation" begin

        include("R_sub_rho.jl")
        ΔS, ΔΘ = diff(salinity[1:2]), diff(temperature[1:2])
        Sₘ, Θₘ = (sum(salinity[1:2]) / 2), (sum(temperature[1:2]) / 2)
        β = haline_contraction(Θₘ, Sₘ, 0, model.buoyancy.model.equation_of_state)
        α = -thermal_expansion(Θₘ, Sₘ, 0, model.buoyancy.model.equation_of_state)

        println(α, β)
        println(sdns.initial_conditions.R_ρ)
        println((β * ΔS) / (α * ΔΘ))
        @test all(sdns.initial_conditions.R_ρ .≈ (β * ΔS) / (α * ΔΘ))

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

end
