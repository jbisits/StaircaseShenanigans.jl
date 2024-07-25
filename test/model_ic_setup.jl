## Random initial conditions
random_no_steps = 3:10
number_of_steps = rand(random_no_steps)
depth_of_steps = reverse(Array(range(domain_extent.Lz+0.2, -0.1, length = number_of_steps)))
salinity = Array(range(34.57, 34.7, length = number_of_steps+1))
temperature = Array(range(-1.5, -1.0, length = number_of_steps+1))

## Setup the model
model = DNSModel(architecture, domain_extent, resolution, diffusivities)

step_ics = StepInitialConditions(model, number_of_steps, depth_of_steps, salinity, temperature)

sdns = StaircaseDNS(model, step_ics)

set_staircase_initial_conditions!(sdns)
