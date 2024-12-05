## Random initial conditions, staircase
random_no_interfaces = 3:10
number_of_interfaces = rand(random_no_interfaces)
depth_of_interfaces = reverse(Array(range(domain_extent.Lz+0.2, -0.1, length = number_of_interfaces)))
salinity_staircase = Array(range(34.57, 34.7, length = number_of_interfaces+1))
temperature_staircase = Array(range(-1.5, -1.0, length = number_of_interfaces+1))

model = DNSModel(architecture, diffusivities, domain_extent, domain_topology, resolution)

step_ics = STStaircaseInitialConditions(model, number_of_interfaces, depth_of_interfaces,
                                        salinity_staircase, temperature_staircase)

sdns_staircase = StaircaseDNS(model, step_ics)

set_initial_conditions!(sdns_staircase)

## Random initial conditions, interface
z = range(domain_extent.Lz, -0.0, length = 20)
depth_of_interface = rand(z[5:15]) # ensures interface is not at extremity of boundary
salinity = Array(range(34.57, 34.7, length = 2))
temperature = Array(range(-1.5, -1.0, length = 2))

model = DNSModel(architecture, diffusivities, domain_extent, domain_topology, resolution)

singel_interface_ics = STSingleInterfaceInitialConditions(model, depth_of_interface, salinity, temperature)

sdns_interface = StaircaseDNS(model, singel_interface_ics)

set_initial_conditions!(sdns_interface)
