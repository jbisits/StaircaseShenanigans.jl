## Simple setup
architecture = CPU()
diffusivities = (ν = 1e-4, κ = (S = 1e-7, T = 1e-5))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = -1.0)
resolution = (Nx = 5, Ny = 5, Nz = 100)

## For testing R_ρ calcualtion
S, Θ = 34.6, -0.5
eos = CustomLinearEquationOfState(Θ, S)
