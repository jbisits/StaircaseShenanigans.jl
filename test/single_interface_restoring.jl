architecture = CPU() # or GPU()
diffusivities = (ν = 1e-5, κ = (S = 1e-8, T = 1e-6))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = -1.0)
resolution = (Nx = 5, Ny = 5, Nz = 70)

number_of_interfaces = 1
depth_of_interfaces = [-0.5]
salinity = [34.58, 34.70]
temperature = [-1.5, 0.5]

model = DNSModel(architecture, domain_extent, resolution, diffusivities)

step_ics = StepInitialConditions(model, number_of_steps, depth_of_steps, salinity, temperature)

sdns = StaircaseDNS(model, step_ics)

set_staircase_initial_conditions!(sdns)

Δx, Δy, Δz = xspacings(model.grid, Center()), yspacings(model.grid, Center()),
                zspacings(model.grid, Center())

ΔV = Δx * Δy * Δz

initial_upper_T_content = sum(interior(sdns.model.tracers.T, :, :, upper)) * ΔV
initial_lower_T_content = sum(interior(sdns.model.tracers.T, :, :, lower)) * ΔV

Δt = 1e-1
stop_time = 50 * 60 # seconds
save_schedule = 10  # seconds
output_path = joinpath(@__DIR__, "output")
simulation = SDNS_simulation_setup(sdns, Δt, stop_time, save_schedule, save_computed_output!,
                                    StaircaseShenanigans.no_velocities!,
                                    S_and_T_tracer_restoring_callbacks!; output_path, max_Δt = 5,
                                    )

run!(simulation)

post_upper_T_content = sum(interior(simulation.model.tracers.T, :, :, upper)) * ΔV
post_lower_T_content = sum(interior(simulation.model.tracers.T, :, :, lower)) * ΔV
