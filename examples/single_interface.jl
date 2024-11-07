using StaircaseShenanigans

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-5, κ = (S = 1e-7, T = 1e-5))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = -1.0)
resolution = (Nx = 5, Ny = 5, Nz = 50)
eos = CustomLinearEquationOfState(-0.5, 34.6)
model_setup = (;architecture, diffusivities, domain_extent, resolution, eos)

## Initial conditions
depth_of_interface = -0.5
salinity = [34.58, 34.70]
temperature = [-1.5, 0.5]
interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature, maintain_interface = true)
velocity_noise = VelocityNoise(0.0, 0.0, 1e-7)

## setup model
sdns = StaircaseDNS(model_setup, interface_ics, velocity_noise)
set_staircase_initial_conditions!(sdns)

## Build simulation
Δt = 1e-1
stop_time = 50 * 60 # seconds
save_schedule = 10  # seconds
output_path = joinpath(@__DIR__, "output")
simulation = SDNS_simulation_setup(sdns, Δt, stop_time, save_schedule, save_computed_output!,
                                    StaircaseShenanigans.save_vertical_velocities!; output_path, max_Δt = 5,
                                    )
## Run
run!(simulation)
