module StaircaseShenanigans

using Oceananigans, Reexport
using SeawaterPolynomials.TEOS10
using SeawaterPolynomials.SecondOrderSeawaterPolynomials
using SeawaterPolynomials: BoussinesqEquationOfState

@reexport using Oceananigans, SeawaterPolynomials.TEOS10,
                SeawaterPolynomials.SecondOrderSeawaterPolynomials

abstract type AbstractStaircaseModel end
abstract type AbstractStaircaseInitialConditions end
abstract type AbstractNoise end

export StaircaseDNS, DNSModel, SDNS_simulation_setup

export StepInitialConditions, SmoothStepInitialConditions, set_staircase_initial_conditions!

include("staircase_model.jl")
include("staircase_initial_conditions.jl")
include("set_staircase_initial_conditions.jl")

end
