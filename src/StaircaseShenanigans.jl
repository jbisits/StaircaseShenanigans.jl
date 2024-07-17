module StaircaseShenanigans

using Oceananigans, Reexport
using SeawaterPolynomials.TEOS10
using SeawaterPolynomials.SecondOrderSeawaterPolynomials
using SeawaterPolynomials: BoussinesqEquationOfState
using TwoLayerDirectNumericalShenanigans: save_tracers!, save_velocities!, no_custom_output!,
                                          simulation_progress, checkpointer_setup!
using TwoLayerDirectNumericalShenanigans: AbstractNoise, SalinityNoise, TemperatureNoise, perturb_tracer

@reexport using Oceananigans, SeawaterPolynomials.TEOS10,
                SeawaterPolynomials.SecondOrderSeawaterPolynomials

abstract type AbstractStaircaseModel end
abstract type AbstractStaircaseInitialConditions end

export StaircaseDNS, DNSModel, SDNS_simulation_setup

export StepInitialConditions, SmoothStepInitialConditions, set_staircase_initial_conditions!

include("staircase_model.jl")
include("staircase_initial_conditions.jl")
include("set_staircase_initial_conditions.jl")

end
