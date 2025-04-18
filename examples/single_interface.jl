using StaircaseShenanigans

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-5, κ = (S = 1e-8, T = 1e-6))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = (-100.0, -90.0))
domain_topology = (x = Periodic, y = Periodic, z = Bounded)
resolution = (Nx = 5, Ny = 5, Nz = 50)
eos = CustomLinearEquationOfState(-0.5, 34.6)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)

## Initial conditions
depth_of_interface = -0.5
salinity = [34.54, 34.70]
temperature = [-1.5, 0.5]
interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature)
velocity_noise = NoiseAtDepth([-0.6, -0.4], VelocityNoise())
tracer_noise = NoiseAtDepth([-0.6, -0.4], TracerNoise())

## setup model
sdns = StaircaseDNS(model_setup, interface_ics, tracer_noise)

# Build simulation
Δt = 1e-1
stop_time = 2 * 60 * 60 # seconds
save_schedule = 60  # seconds
output_path = joinpath(@__DIR__, "output_non_periodic")
simulation = SDNS_simulation_setup(sdns, stop_time, save_computed_output!, save_vertical_velocities!;
                                    Δt, save_schedule,
                                    output_path, max_Δt = 5)
## Run
run!(simulation)
compute_R_ρ!(simulation.output_writers[:computed_output].filepath,
             simulation.output_writers[:tracers].filepath, eos)
