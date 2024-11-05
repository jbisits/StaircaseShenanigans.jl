model = DNSModel(architecture, domain_extent, resolution, diffusivities, eos)

## Initial conditions
number_of_interfaces = 4
depth_of_interfaces = [-0.2, -0.4, -0.6, -0.8]
salinity = [34.57, 34.6, 34.63, 34.66, 34.69]
temperature = [-1.5, -1.0, -0.5, 0.0, 0.5]

step_ics = STStaircaseInitialConditions(model, number_of_interfaces, depth_of_interfaces, salinity, temperature)
interface_ics = STSingleInterfaceInitialConditions(model, depth_of_interfaces[2], salinity[1:2], temperature[1:2])
