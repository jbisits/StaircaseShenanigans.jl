module StaircaseShenanigans

using Oceananigans, Reexport, Printf
using Oceananigans: seawater_density
using Oceananigans.BoundaryConditions: update_boundary_condition!
using Oceananigans.Operators: ℑzᵃᵃᶠ
using Oceananigans.Fields: condition_operand
import Oceananigans.BackgroundField
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
abstract type BackgroundFunction end

export StaircaseDNS, PeriodicStaircaseDNS, DNSModel, SDNS_simulation_setup

export STStaircaseInitialConditions, StaircaseICs, SmoothSTStaircaseInitialConditions,
       STSingleInterfaceInitialConditions, SingleInterfaceICs,
       PeriodicSTSingleInterfaceInitialConditions, PeriodoicSingleInterfaceICs,
       set_staircase_initial_conditions!

export BackgroundTanh, BackgroundLinear, tanh_background, linear_background

export AbstractNoise, VelocityNoise, TracerNoise

export OuterStairMask, OuterStairTargets, OuterMask, OuterTargets, ExponentialTarget

export CustomLinearRoquetSeawaterPolynomial, CustomLinearEquationOfState

export save_computed_output!, save_all_velocities!, save_vertical_velocities!

export compute_R_ρ!

export animate_tracers, animate_density, visualise_initial_conditions, visualise_initial_density,
       animate_tracers_anomaly, animate_density_anomaly

include("staircase_initial_conditions.jl")
include("staircase_noise.jl")
include("staircase_diagnostics.jl")
include("staircase_model.jl")
include("set_staircase_initial_conditions.jl")
include("staircase_background.jl")
include("staircase_restoring.jl")
include("alt_lineareos.jl")
include("makie_plotting_functions.jl")

end
