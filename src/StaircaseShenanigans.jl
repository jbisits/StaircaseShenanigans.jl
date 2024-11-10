module StaircaseShenanigans

using Oceananigans, Reexport, Printf
using Oceananigans: seawater_density
using Oceananigans.BoundaryConditions: update_boundary_condition!
using Oceananigans.Operators: ℑzᵃᵃᶠ
using SeawaterPolynomials.TEOS10
using SeawaterPolynomials.SecondOrderSeawaterPolynomials
using SeawaterPolynomials: BoussinesqEquationOfState
using SeawaterPolynomials: thermal_expansion, haline_contraction, ρ
using GibbsSeaWater: gsw_alpha, gsw_beta
using NCDatasets, JLD2
using Statistics

import SeawaterPolynomials.SecondOrderSeawaterPolynomials: RoquetSeawaterPolynomial

@reexport using Oceananigans, SeawaterPolynomials.TEOS10,
                SeawaterPolynomials.SecondOrderSeawaterPolynomials

abstract type AbstractStaircaseModel end
abstract type AbstractInitialConditions end
abstract type AbstractNoise end

export StaircaseDNS, PeriodicStaircaseDNS, DNSModel, SDNS_simulation_setup

export STStaircaseInitialConditions, StaircaseICs, SmoothSTStaircaseInitialConditions,
       STSingleInterfaceInitialConditions, SingleInterfaceICs,
       PeriodicSTSingleInterfaceInitialConditions, PeriodoicSingleInterfaceICs,
       set_staircase_initial_conditions!

export tanh_background

export AbstractNoise, VelocityNoise, TracerNoise

export OuterStairMask, OuterStairTargets, OuterMask, OuterTargets, ExponentialTarget

export CustomLinearRoquetSeawaterPolynomial, CustomLinearEquationOfState

export save_computed_output!

export animate_tracers, animate_density, visualise_initial_conditions, visualise_initial_density

include("staircase_initial_conditions.jl")
include("staircase_noise.jl")
include("staircase_model.jl")
include("set_staircase_initial_conditions.jl")
include("staircase_restoring.jl")
include("alt_lineareos.jl")
include("makie_plotting_functions.jl")

end
