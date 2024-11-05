using StaircaseShenanigans
using StaircaseShenanigans: Oceananigans.Operators.ℑzᵃᵃᶠ

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-5, κ = (S = 1e-7, T = 1e-5))
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

# @inline tracer_flux_change(i, j, grid, clock, model_fields, ΔC) =
#     @inbounds - ΔC * (model_fields.w[i, j, grid.Nz] * ℑzᵃᵃᶠ(model_fields.T[i, j, grid.Nz])) /
#                     volume(i, j, k, grid)

@inline tracer_flux_change(i, j, grid, clock, model_fields, ΔT) =
    @inbounds + ΔT + model_fields.T[i, j, grid.Nz] #ℑzᵃᵃᶠ(i, j, grid.Nz, model_fields.T)


@inline w_top_bottom(i, j, grid, clock, model_fields) =
    @inbounds model_fields.w[i, j, 1]
@inline w_bottom_top(i, j, grid, clock, model_fields) =
    @inbounds model_fields.w[i, j, grid.Nz+1]

ΔT = 2
T_reentry = ValueBoundaryCondition(tracer_flux_change, discrete_form=true, parameters = ΔT)
T_bcs = FieldBoundaryConditions(bottom = T_reentry)
w_top = OpenBoundaryCondition(w_top_bottom, discrete_form=true)
w_bottom = OpenBoundaryCondition(w_bottom_top, discrete_form=true)
# w_top = OpenBoundaryCondition(1e-5)
# w_bottom = OpenBoundaryCondition(1e-5)
w_bcs = FieldBoundaryConditions(top = w_top,  bottom = w_bottom)
boundary_conditions = (T=T_bcs, w= w_bcs)
## model
model = DNSModel(architecture, domain_extent, resolution, diffusivities, eos; boundary_conditions)

## Set initial conditions
step_ics = StepInitialConditions(model, number_of_interfaces, depth_of_interfaces, salinity, temperature)

sdns = StaircaseDNS(model, step_ics)

set_staircase_initial_conditions!(sdns)

w_noise = randn(size(sdns.model.velocities.w)) * 1e-7
set!(sdns.model, w = w_noise)

## Build simulation
Δt = 1e-1
stop_time = 50 * 60 # seconds
save_schedule = 10  # seconds
output_path = joinpath(@__DIR__, "output")
simulation = SDNS_simulation_setup(sdns, Δt, stop_time, save_schedule, save_computed_output!, StaircaseShenanigans.save_vertical_velocities!; output_path, max_Δt = 5,
                                    )

# ##
# function modify_halo!(simulation, parameters)

#     T = simulation.model.tracers.T
#     T[1:5, 1:5, -2:1] .+= parameters.ΔT
#     T[1:5, 1:5, 50:53] .-= parameters.ΔT
#     return nothing
# end
# simulation.callbacks[:T_bc_restore] = Callback(modify_halo!, parameters = (ΔT = 2,))
## Run
run!(simulation)
