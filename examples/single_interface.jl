using StaircaseShenanigans

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-4, κ = (S = 1e-7, T = 1e-5))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = -1.0)
resolution = (Nx = 5, Ny = 5, Nz = 50)

## Initial conditions
number_of_interfaces = 1
depth_of_interfaces = [-0.5]
salinity = [34.58, 34.70]
temperature = [-1.5, 0.5]

## Setup the model
eos = CustomLinearEquationOfState(0, 34.6)

## restoring
rate = 1/10
mask = OuterMask(-0.25, -0.75)
S_target = OuterTargets(34.58, 34.7, mask)
T_target = OuterTargets(-1.5, 0.5, mask)
T_restoring = Relaxation(; rate, mask, target = T_target)
S_restoring = Relaxation(; rate, mask, target = S_target)
forcing = (S = S_restoring, T = T_restoring)

## model
model = DNSModel(architecture, domain_extent, resolution, diffusivities, eos)

## Set initial conditions
step_ics = StepInitialConditions(model, number_of_interfaces, depth_of_interfaces, salinity, temperature)

sdns = StaircaseDNS(model, step_ics)

set_staircase_initial_conditions!(sdns)

## Build simulation
Δt = 1e-1
stop_time = 50 * 60 # seconds
save_schedule = 10  # seconds
output_path = joinpath(@__DIR__, "output")
simulation = SDNS_simulation_setup(sdns, Δt, stop_time, save_schedule, save_computed_output!,
                                    StaircaseShenanigans.no_velocities!,
                                    S_and_T_tracer_restoring_callbacks!; output_path, max_Δt = 5,
                                    )

## Run
run!(simulation)
