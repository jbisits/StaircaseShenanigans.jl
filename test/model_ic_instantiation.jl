## single interface
depth_of_interface =  -0.6
salinity = [34.57, 34.69]
temperature = [-1.5, 0.5]

## staircase
number_of_interfaces = 2
depth_of_interfaces =  [-0.4, -0.6]
salinities = [34.57, 34.69, 34.8]
temperatures = [-1.5, 0.5, 1.0]

smoothing = (NoSmoothing, TanhInterfaceThickness(), TanhInterfaceThickness(0.05), TanhInterfaceThickness(0.01, 0.05))
background = (NoBackground, BackgroundTanh(), BackgroundLinear(), BackgroundStep())
