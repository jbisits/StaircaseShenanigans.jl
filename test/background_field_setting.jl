
model = DNSModel(architecture, diffusivities, domain_extent, domain_topology, resolution)

singel_interface_ics = STSingleInterfaceInitialConditions(model, depth_of_interface, salinity, temperature,
                                                            background_state = BackgroundLinear())

sdns_interface_linear_bf = StaircaseDNS(model, singel_interface_ics)

set_initial_conditions!(sdns_interface_linear_bf)

singel_interface_ics = STSingleInterfaceInitialConditions(model, depth_of_interface, salinity, temperature,
                                                            background_state = BackgroundTanh())

sdns_interface_tanh_bf = StaircaseDNS(model, singel_interface_ics)
