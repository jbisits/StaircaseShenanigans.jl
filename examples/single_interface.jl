using StaircaseShenanigans

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-5, κ = (S = 1e-8, T = 1e-6))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = -1.0)
resolution = (Nx = 5, Ny = 5, Nz = 50)

## Initial conditions
number_of_interfaces = 1
depth_of_interfaces = [-0.5]
salinity = [34.58, 34.70]
temperature = [-1.5, 0.5]

## Setup the model
eos = CustomLinearEquationOfState(-0.5, 34.6)

## restoring
# rate = 1/10
# mask = OuterMask(-0.25, -0.75)
# S_target = OuterTargets(34.58, 34.7, mask)
# T_target = OuterTargets(-1.5, 0.5, mask)
# T_restoring = Relaxation(; rate, mask, target = T_target)
# S_restoring = Relaxation(; rate, mask, target = S_target)
# forcing = (S = S_restoring, T = T_restoring)

# discrete forcing
T= model.tracers.T
upper_region(i, j, k, grid, T) = k < grid.Nz * 0.4
lower_region(i, j, k, grid, T) = k > grid.Nz * 0.6
T_upper = condition_operand(T, upper_region, 0)
∫T_upper = Field(Integral(T_upper))
compute!(∫T_upper)
T_lower = condition_operand(T, lower_region, 0)
∫T_lower = Field(Integral(T_upper))
compute!(∫T_lower)
initial_T_content = (upper = ∫T_upper.data[1, 1, 1], lower = ∫T_lower.data[1, 1, 1])

function T_forcing_func(i, j, k, grid, clock, model_fields, initial_T_content)

    T = model_fields.T
    upper_region(i, j, k, grid, T) = k < grid.Nz * 0.4
    lower_region(i, j, k, grid, T) = k > grid.Nz * 0.6
    T_upper = condition_operand(T, upper_region, 0)
    ∫T_upper = Field(Integral(T_upper))
    ΔT_upper = initial_T_content.upper[1, 1, 1] - ∫T_upper[1, 1, 1]
    T_lower = condition_operand(T, lower_region, 0)
    ∫T_lower = Field(Integral(T_lower))
    ΔT_lower = initial_T_content.lower[1, 1, 1] - ∫T_lower[1, 1, 1]
    return k < grid.Nz * 0.1 ? 0.1 * ΔT_upper[1, 1, 1] / 1024 : k > grid.Nz * 0.9 ?
                                                                0.1 * ΔT_lower[1, 1, 1] / 1024 : 0
end
T_forcing = Forcing(T_forcing_func, discrete_form = true, parameters = initial_T_content)
forcing = (T = T_forcing, )

## model
model = DNSModel(architecture, domain_extent, resolution, diffusivities, eos; forcing)

## Set initial conditions
step_ics = StepInitialConditions(model, number_of_interfaces, depth_of_interfaces, salinity, temperature)

sdns = StaircaseDNS(model, step_ics)

set_staircase_initial_conditions!(sdns)

## Build simulation
Δt = 1e-1
stop_time = 50 * 60 # seconds
save_schedule = 10  # seconds
output_path = joinpath(@__DIR__, "output")
simulation = SDNS_simulation_setup(sdns, Δt, stop_time, save_schedule, save_computed_output!; output_path, max_Δt = 5,
                                    )

## Run
run!(simulation)
