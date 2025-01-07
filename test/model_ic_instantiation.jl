## single interface
depth_of_interface =  -0.6
salinity = [34.57, 34.69]
temperature = [-1.5, 0.5]

smoothing = (NoSmoothing, TanhInterfaceSteepness(), TanhInterfaceSteepness(100.0), TanhInterfaceSteepness(1000.0, 200.0))
background = (NoBackground, BackgroundTanh(), BackgroundLinear(), BackgroundStep())
