using StaircaseShenanigans, GibbsSeaWater

restart = true

architecture = CPU()
Pr = 7   # Prandtl
τ = 0.05 # diff ratio
ν = 2.5e-6 # set this get the others
diffusivities = diffusivities_from_ν(ν; τ, Pr)
domain_extent = (Lx=0.05, Ly=0.05, Lz=-1.0)
domain_topology = (x = Periodic, y = Periodic, z = Bounded)
resolution = (Nx=5, Ny=5, Nz=100)
ρ₀ = gsw_rho(34.7, 0.5, 0)
eos = TEOS10EquationOfState(reference_density = ρ₀)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)
model = DNSModel(model_setup...; TD = VerticallyImplicitTimeDiscretization())

number_of_interfaces = 3
depth_of_interfaces = [-0.25, -0.5, -0.75]
#### nonlinear eos
salinity = [34.58, 34.61, 34.65, 34.7]
temperature = [-1.51, -0.87, -0.19, 0.51]
staircase_ics = StaircaseICs(model, number_of_interfaces, depth_of_interfaces, salinity, temperature)

initial_noise = (velocities = VelocityNoise(1e-2), tracers = TracerNoise(1e-4, 1e-2))

## setup model
sdns = StaircaseDNS(model, staircase_ics; initial_noise)

## Build simulation
stop_time = Int(1 * 60 * 60) # seconds
initial_state = "step" # can update if smoothing is added
output_path = joinpath(@__DIR__, "rundown_$(round(staircase_ics.R_ρ[2], digits = 2))", initial_state)
save_schedule = 60
checkpointer_time_interval = 60 * 60 # seconds
Δt = 1e-1
simulation = SDNS_simulation_setup(sdns, stop_time, save_computed_output!,
                                   save_vertical_velocities!;
                                   output_path,
                                   save_schedule,
                                   checkpointer_time_interval,
                                   overwrite_saved_output = restart,
                                   Δt)
## Run
restart ? nothing : simulation.stop_time = Int(10 * 60 * 60) # update to pickup from a checkpoint
pickup = restart ? false : readdir(simulation.output_writers[:checkpointer].dir, join = true)[1]
run!(simulation; pickup)
