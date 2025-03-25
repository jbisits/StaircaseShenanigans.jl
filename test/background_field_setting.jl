depth_of_interface =  -0.6
salinity = [34.57, 34.7]
temperature = [-1.5, 0.5]

model = DNSModel(architecture, diffusivities, domain_extent, domain_topology, resolution)

singel_interface_ics = STSingleInterfaceInitialConditions(model_setup.eos, depth_of_interface, salinity, temperature,
                                                            background_state = BackgroundLinear())

sdns_interface_linear_bf = StaircaseDNS(model_setup, singel_interface_ics, nothing)

set_initial_conditions!(sdns_interface_linear_bf)

singel_interface_ics = STSingleInterfaceInitialConditions(model_setup.eos, depth_of_interface, salinity, temperature,
                                                            background_state = BackgroundTanh())

sdns_interface_tanh_bf = StaircaseDNS(model_setup, singel_interface_ics, nothing)

set_initial_conditions!(sdns_interface_tanh_bf)
