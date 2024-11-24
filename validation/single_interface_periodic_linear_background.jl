using StaircaseShenanigans

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-5, κ = (S = 1e-7, T = 1e-5))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = -1.0)
resolution = (Nx = 5, Ny = 5, Nz = 50)
eos = CustomLinearEquationOfState(-0.5, 34.6)
model_setup = (;architecture, diffusivities, domain_extent, resolution, eos)

## Initial conditions
depth_of_interface = -0.5
salinity = [34.56, 34.70]
temperature = [-1.5, 0.5]

#### Background linear state with tanh smoothing
interface_ics = PeriodoicSingleInterfaceICs(eos, depth_of_interface, salinity,
                                            background_state = BackgroundLinear(),
                                            interface_smoothing = Tanh)

## setup model
sdns = StaircaseDNS(model_setup, interface_ics, nothing)

## Build simulation
Δt = 1e-1
stop_time = 4 * 60 * 60 # seconds
save_schedule = 30  # seconds
output_path = joinpath(@__DIR__, "output_periodic_tanh_interface_linear_background")
simulation = SDNS_simulation_setup(sdns, stop_time, save_computed_output!,
                                    save_vertical_velocities!;
                                    Δt, save_schedule,
                                    output_path, max_Δt = 5)
## Run
run!(simulation)
compute_R_ρ!(simulation.output_writers[:computed_output].filepath,
             simulation.output_writers[:tracers].filepath, eos)
