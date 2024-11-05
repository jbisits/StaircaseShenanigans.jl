using StaircaseShenanigans

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-5, κ = (S = 1e-7, T = 1e-5))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = -1.0)
resolution = (Nx = 5, Ny = 5, Nz = 50)

## Initial conditions
number_of_interfaces = 1
depth_of_interfaces = [-0.5]
salinity = [34.58, 34.70]
temperature = [-1.5, 0.5]

## Setup the model
eos = CustomLinearEquationOfState(-0.5, 34.6)


## model
model = DNSModel(architecture, domain_extent, resolution, diffusivities, eos)

## Set initial conditions
step_ics = StepInitialConditions(model, number_of_interfaces, depth_of_interfaces, salinity, temperature)

sdns = StaircaseDNS(model, step_ics)

set_staircase_initial_conditions!(sdns)

w_noise = randn(size(sdns.model.velocities.w)) * 1e-7
set!(sdns.model, w = w_noise)

## Build simulation
Δt = 1e-1
stop_time = 50 * 60 # seconds
save_schedule = 10  # seconds
output_path = joinpath(@__DIR__, "output/triple_periodic")
simulation = SDNS_simulation_setup(sdns, Δt, stop_time, save_schedule, save_computed_output!, StaircaseShenanigans.save_vertical_velocities!; output_path, max_Δt = 5,
                                    )

# ##
# function modify_halo!(simulation, parameters)

#     T = simulation.model.tracers.T
#     T[1:5, 1:5, -2:1] .+= parameters.ΔT
#     T[1:5, 1:5, 50:53] .-= parameters.ΔT
#     return nothing
# end
# simulation.callbacks[:T_bc_restore] = Callback(modify_halo!, parameters = (ΔT = 2,))
## Run
run!(simulation)
