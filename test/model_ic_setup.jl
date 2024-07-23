architecture = CPU()
diffusivities = (ν = 1e-4, κ = (S = 1e-7, T = 1e-5))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = -1.0)
resolution = (Nx = 5, Ny = 5, Nz = 100)

## Setup the model
model = DNSModel(architecture, domain_extent, resolution, diffusivities)

## Initial conditions
random_no_steps = 3:10
number_of_steps = rand(random_no_steps)
depth_of_steps = reverse(Array(range(domain_extent.Lz+0.2, -0.1, length = number_of_steps)))
salinity = Array(range(34.57, 34.7, length = number_of_steps+1))
temperature = Array(range(-1.5, -1.0, length = number_of_steps+1))

step_ics = StepInitialConditions(number_of_steps, depth_of_steps, salinity, temperature)

sdns = StaircaseDNS(model, step_ics)

set_staircase_initial_conditions!(sdns)
