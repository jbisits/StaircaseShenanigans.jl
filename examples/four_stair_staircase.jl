using StaircaseShenanigans

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-4, κ = (S = 1e-7, T = 1e-5))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = -1.0)
resolution = (Nx = 5, Ny = 5, Nz = 100)

## Initial conditions
number_of_steps = 4
depth_of_steps = [-0.2, -0.4, -0.6, -0.8]
salinity = [34.57, 34.60, 34.63, 34.66, 34.69]
temperature = [-1.5, -1.0, -0.5, 0.0, 0.5]

## Setup the model
rate = 1/40
# TODO: make this a nice function that can take different stairs and restore
# How to restore to seperate values in each stair?
test_mask(x, y, z) = z > -0.2 ? 1 : 0
function another_mask(x, y, z)
    if z > -0.2
        1
    elseif z < -0.8
        (1 + 5 / 33)
    else
        0
    end
end
target_salinity = 33
S_relaxation = Relaxation(; rate, mask = another_mask, target = target_salinity)
relaxation = (S = S_relaxation,)
model = DNSModel(architecture, domain_extent, resolution, diffusivities; relaxation)

## Set initial conditions
step_ics = StepInitialConditions(model, number_of_steps, depth_of_steps, salinity, temperature)

sdns = StaircaseDNS(model, step_ics)

set_staircase_initial_conditions!(sdns)

## Build simulation
Δt = 1e-1
stop_time = 2 * 60 # seconds
save_schedule = 10  # seconds
output_path = joinpath(@__DIR__, "output")
simulation = SDNS_simulation_setup(sdns, Δt, stop_time, save_schedule, save_computed_output!; output_path)

## Run
run!(simulation)
