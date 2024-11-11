"Condition for lower quarter of domain."
lower_quarter(i, j, k, grid, c) = k < grid.Nz / 4
"Condition for lower quarter of domain."
upper_quarter(i, j, k, grid, c) = 3 * grid.Nz / 4 < k < grid.Nz

"""
    function compute_R_ρ!(computed_output::AbstractString, tracers::AbstractString, eos)
From saved `tracers` output get the averaged salinity and temperature and lower and upper
quarter of the domain, compute the density ratio and save to `computed_output`.
**Note:** this method *needs the equation of state* used in the model so best practice is to
compute this at the end of a script so everything can be easily accessed (see a singlge
interface example).
"""
function compute_R_ρ!(computed_output::AbstractString, tracers::AbstractString, eos)

    ds = NCDataset(tracers)

    S_u = S_g = ds[:Sᵤ_mean][:]
    S_l = S_f = ds[:Sₗ_mean][:]
    T_u = T_f = ds[:Tᵤ_mean][:]
    T_l = T_g = ds[:Tₗ_mean][:]

    eos_vec = fill(eos, length(S_u))
    ρ_u = ρ.(T_u, S_u, 0, eos_vec)
    ρ_l = ρ.(T_l, S_l, 0, eos_vec)
    ρ_f = ρ.(T_f, S_f, 0, eos_vec)
    ρ_g = ρ.(T_g, S_g, 0, eos_vec)

    R_ρ = @. (0.5 * (ρ_f - ρ_u) + 0.5 * (ρ_l - ρ_g)) / (0.5 * (ρ_f - ρ_l) + 0.5 * (ρ_u - ρ_g))

    NCDataset(computed_output, "a") do ds2
        defVar(ds2, "R_ρ", R_ρ, ("time",),
                attrib = Dict("long_name" => "Density ratio"))
    end

    close(ds)

    return nothing
end
