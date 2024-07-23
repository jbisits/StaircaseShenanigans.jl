using StaircaseShenanigans
using Test

@testset "StaircaseShenanigans.jl" begin

    @testset begin
        include("model_ic_setup.jl")

        @test all(unique(interior(sdns.model.tracers.S, 1, 1, :)) .== reverse(salinity))
        @test all(unique(interior(sdns.model.tracers.T, 1, 1, :)) .== reverse(temperature))

    end
end
