model = DNSModel(architecture, domain_extent, resolution, diffusivities, eos)

## Initial conditions
number_of_steps = 4
depth_of_steps = [-0.2, -0.4, -0.6, -0.8]
salinity = [34.57, 34.6, 34.63, 34.66, 34.69]
temperature = [-1.5, -1.0, -0.5, 0.0, 0.5]

step_ics = StepInitialConditions(model, number_of_steps, depth_of_steps, salinity, temperature)

sdns = StaircaseDNS(model, step_ics)

set_staircase_initial_conditions!(sdns)
