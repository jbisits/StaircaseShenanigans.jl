module StaircaseShenanigans

using Oceananigans, Reexport, Printf
using Oceananigans: seawater_density
using SeawaterPolynomials.TEOS10
using SeawaterPolynomials.SecondOrderSeawaterPolynomials
using SeawaterPolynomials: BoussinesqEquationOfState
using NCDatasets, JLD2

@reexport using Oceananigans, SeawaterPolynomials.TEOS10,
                SeawaterPolynomials.SecondOrderSeawaterPolynomials

abstract type AbstractStaircaseModel end
abstract type AbstractStaircaseInitialConditions end
abstract type AbstractNoise end

export StaircaseDNS, DNSModel, SDNS_simulation_setup

export StepInitialConditions, SmoothStepInitialConditions, set_staircase_initial_conditions!

export save_computed_output!

include("staircase_model.jl")
include("staircase_initial_conditions.jl")
include("set_staircase_initial_conditions.jl")

end
