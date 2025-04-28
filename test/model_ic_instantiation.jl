## single interface
depth_of_interface =  -0.6
salinity = [34.57, 34.69]
temperature = [-1.5, 0.5]

## staircase
depth_of_interfaces =  -0.6
salinities = [34.57, 34.69]
temperatures = [-1.5, 0.5]

smoothing = (NoSmoothing, TanhInterfaceThickness(), TanhInterfaceThickness(0.05), TanhInterfaceThickness(0.01, 0.05))
background = (NoBackground, BackgroundTanh(), BackgroundLinear(), BackgroundStep())
