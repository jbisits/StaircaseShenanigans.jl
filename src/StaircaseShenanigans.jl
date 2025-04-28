module StaircaseShenanigans

using Oceananigans, Reexport, Printf
using Oceananigans: seawater_density
using Oceananigans.BuoyancyFormulations: buoyancy_frequency
using Oceananigans.Fields: condition_operand
import Oceananigans.BackgroundField
using Oceananigans.TurbulenceClosures: cell_diffusion_timescale, ExplicitTimeDiscretization
using SeawaterPolynomials.TEOS10
using SeawaterPolynomials.SecondOrderSeawaterPolynomials
using SeawaterPolynomials: BoussinesqEquationOfState
using SeawaterPolynomials: thermal_expansion, haline_contraction, total_density
using GibbsSeaWater: gsw_alpha, gsw_beta
using NCDatasets, JLD2
using Statistics
using Oceanostics: KineticEnergyDissipationRate, KineticEnergy, PotentialEnergy, BuoyancyProductionTerm

import SeawaterPolynomials.SecondOrderSeawaterPolynomials: RoquetSeawaterPolynomial

@reexport using Oceananigans, SeawaterPolynomials.TEOS10,
                SeawaterPolynomials.SecondOrderSeawaterPolynomials

abstract type AbstractStaircaseModel end
abstract type AbstractInitialConditions end
abstract type AbstractNoise end
abstract type AbstractBackgroundFunction end
abstract type AbstractInterfaceSmoothing end

export StaircaseDNS, DNSModel, SDNS_simulation_setup, enhance_κₛ, enhance_κₜ,
       diffusivities_from_ν, diffusivities_from_κₜ

export STStaircaseInitialConditions, StaircaseICs,
       STSingleInterfaceInitialConditions, SingleInterfaceICs,
       set_initial_conditions!

export TanhInterfaceSmoothing, Tanh, TanhInterfaceThickness, NoSmoothing

export BackgroundTanh, BackgroundLinear, BackgroundStep, NoBackground,
        tanh_background, linear_background, step_background

export jump_periodic_boundary_conditions

export AbstractNoise, VelocityNoise, TracerNoise, NoiseAtDepth

export OuterStairMask, OuterStairTargets, OuterMask, OuterTargets, ExponentialTarget

export CustomLinearRoquetSeawaterPolynomial, CustomLinearEquationOfState

export save_computed_output!, save_all_velocities!, save_vertical_velocities!

export compute_R_ρ!, save_diagnostics!, save_diagnostic!, update_diagnostic!

export animate_tracers, animate_density, visualise_initial_conditions, visualise_initial_density,
       animate_tracers_anomaly, animate_density_anomaly, animate_profile_in_S_Θ_space

include("staircase_initial_conditions.jl")
include("single_interfaces_initial_conditions.jl")
include("interface_smoothing.jl")
include("staircase_background.jl")
include("staircase_jumpperiodic.jl")
include("staircase_noise.jl")
include("staircase_diagnostics.jl")
include("staircase_model.jl")
include("staircase_simulation.jl")
include("set_initial_conditions.jl")
include("staircase_restoring.jl")
include("alt_lineareos.jl")
include("makie_plotting_functions.jl")

end
