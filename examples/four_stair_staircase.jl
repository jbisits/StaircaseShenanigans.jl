using StaircaseShenanigans

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-4, κ = (S = 1e-7, T = 1e-5))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = -1.0)
resolution = (Nx = 5, Ny = 5, Nz = 50)

## Initial conditions
number_of_interfaces = 4
depth_of_interfaces = [-0.2, -0.4, -0.6, -0.8]
salinity = [34.55, 34.60, 34.65, 34.70, 34.75]
temperature = [-1.5, -1.0, -0.5, 0.0, 0.5]

## Setup the model
eos = CustomLinearEquationOfState(0, 34.6)
model = DNSModel(architecture, domain_extent, resolution, diffusivities, eos)

## Set initial conditions
step_ics = StaircaseICs(model, number_of_interfaces, depth_of_interfaces, salinity, temperature)

sdns = StaircaseDNS(model, step_ics)

set_staircase_initial_conditions!(sdns)

## Build simulation
Δt = 1e-1
stop_time = 10 * 60 # seconds
save_schedule = 10  # seconds
output_path = joinpath(@__DIR__, "output")
simulation = SDNS_simulation_setup(sdns, Δt, stop_time, save_schedule, save_computed_output!; output_path)

## Run
run!(simulation)
