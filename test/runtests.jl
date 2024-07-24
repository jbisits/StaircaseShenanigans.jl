using StaircaseShenanigans, GibbsSeaWater
using Test

@testset "StaircaseShenanigans.jl" begin

    @testset "Setting initial conditions" begin
        include("model_ic_setup.jl")

        @test all(unique(interior(sdns.model.tracers.S, 1, 1, :)) .== reverse(salinity))
        @test all(unique(interior(sdns.model.tracers.T, 1, 1, :)) .== reverse(temperature))

    end

    @testset "CustomLienarEquationOfState initialisation" begin
        for _ in 1:5
            Θ = rand(range(-1.5, 20.5, length = 20))
            S = rand(range(34.5, 34.7, length = 20))
            true_coefficients = gsw_alpha(S, Θ, 0), gsw_beta(S, Θ, 0)
            custom_linear = CustomLinearEquationOfState(Θ, S)
            sp = custom_linear.seawater_polynomial
            linear_coefficients = sp.R₀₁₀, sp.R₁₀₀
            other_coefficients = sp.R₁₀₁, sp.R₂₀₀, sp.R₀₂₀, sp.R₁₁₀
            @test sp isa SecondOrderSeawaterPolynomial
            @test all(other_coefficients .== 0)
            @test all(linear_coefficients .== true_coefficients)
        end
    end

end
