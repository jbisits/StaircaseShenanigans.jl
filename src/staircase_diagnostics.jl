"Condition for lower quarter of domain."
lower_quarter(i, j, k, grid, c) = k < 0.25 * grid.Nz
"Condition for lower quarter of domain."
upper_quarter(i, j, k, grid, c) = 0.75 * grid.Nz < k < grid.Nz

"Condition for middle fifth of upper half of domain."
upper_middle_fifth(i, j, k, grid, c) = 0.2 * grid.Nz < k < 0.4 * grid.Nz
"Condition for middle fifth of lower half of domain."
lower_middle_fifth(i, j, k, grid, c) = 0.6 * grid.Nz < k < 0.8 * grid.Nz
"""
    function compute_R_ρ!(computed_output::AbstractString, tracers::AbstractString, eos)
From saved `tracers` output get the averaged salinity and temperature and lower and upper
quarter of the domain, compute the density ratio and save to `computed_output`.
**Note:** this method *needs the equation of state* used in the model so best practice is to
compute this at the end of a script so everything can be easily accessed (see a singlge
interface example).
"""
function compute_R_ρ!(computed_output::AbstractString, tracers::AbstractString,
                      upper::Tuple, lower::Tuple, eos; i = "")

    interface_depth = NCDataset(computed_output) do co
                          co.attrib[:interface_depth]
                      end
    ds = NCDataset(tracers)

    z = ds[:z_aac][:]

    upper_range = findall(upper[1] .< z .< upper[2])
    lower_range = findall(lower[1] .< z .< lower[2])

    S_u = S_g = reshape(mean(ds[:S_ha][upper_range, :], dims = 1), :)
    S_l = S_f = reshape(mean(ds[:S_ha][lower_range, :], dims = 1), :)
    T_u = T_f = reshape(mean(ds[:T_ha][upper_range, :], dims = 1), :)
    T_l = T_g = reshape(mean(ds[:T_ha][lower_range, :], dims = 1), :)

    eos_vec = fill(eos, length(S_u))
    ρ_u = total_density.(T_u, S_u, interface_depth, eos_vec)
    ρ_l = total_density.(T_l, S_l, interface_depth, eos_vec)
    ρ_f = total_density.(T_f, S_f, interface_depth, eos_vec)
    ρ_g = total_density.(T_g, S_g, interface_depth, eos_vec)

    R_ρ = @. (0.5 * (ρ_f - ρ_u) + 0.5 * (ρ_l - ρ_g)) / (0.5 * (ρ_f - ρ_l) + 0.5 * (ρ_u - ρ_g))

    NCDataset(computed_output, "a") do ds2
        if haskey(ds2, "R_ρ"*i)
            # rename so if picking up can saved the whole array
            renameVar(ds2, "R_ρ"*i, "R_ρ_prior_cp"*i)
            defVar(ds2, "R_ρ"*i, R_ρ, ("time",),
                    attrib = Dict("long_name" => "Density ratio at interface "*i))
        else
            defVar(ds2, "R_ρ"*i, R_ρ, ("time",),
                    attrib = Dict("long_name" => "Density ratio at interface "*i))
        end
    end

    close(ds)

    return nothing
end
function compute_R_ρ!(computed_output::AbstractString, tracers::AbstractString,
                      upper::Array{Tuple{Float64, Float64}}, lower::Array{Tuple{Float64, Float64}}, eos)

    for i ∈ eachindex(upper)
        compute_R_ρ!(computed_output, tracers, upper[i], lower[i], eos, i = "$i")
    end

    return nothing
end
"""
    function save_diagnostics!(diagnostics_file::AbstractString, tracers::AbstractString,
                            computed_output::AbstractString; group = nothing)
Save the diagnostics (in this script) to `diagnostics_file`. These diagnostics
are currently only able to be computed from output that is saved in netcdf format (i.e. `.nc` files)
but they returned output is in `.jld2` format. The kwarg `group` is for creating a `group` in
`diagnostics_file`.
"""
function save_diagnostics!(diagnostics_file::AbstractString, tracers::AbstractString,
                           computed_output::AbstractString, velocities::AbstractString;
                           group = nothing,
                           interface_offset = 4)

    group = isnothing(group) ? "" : group[end] == '/' ? group : group * "/" # creates a group in the saved output.

    if isfile(diagnostics_file)

        group_keys = jldopen(diagnostics_file) do f; keys(f); end
        saved_keys = jldopen(diagnostics_file) do f; [keys(f[g]) for g ∈ group_keys]; end
        saved_keys = vcat(saved_keys...)

        if group[end-1] ∉ group_keys

            "dims" ∈ group_keys ? nothing : dimensions!(diagnostics_file, computed_output)
            save_computed_output!(diagnostics_file, computed_output, group)
            potential_and_background_potential_energy!(diagnostics_file, computed_output, tracers, group)
            φ_interface_flux!(diagnostics_file, tracers, :S, group)
            φ_interface_flux!(diagnostics_file, tracers, :T, group)
            ha_φ_flux!(diagnostics_file, tracers, :S_ha, group)
            ha_φ_flux!(diagnostics_file, tracers, :T_ha, group)
            ∫gρw!(diagnostics_file, computed_output, velocities, group)
            initial_non_dim_numbers!(diagnostics_file, computed_output, group)
            save_horizontally_averaged_fields!(diagnostics_file, computed_output, tracers, group)
            compute_Ẽ!(diagnostics_file, computed_output, tracers, group, interface_offset)
            interface_thickness!(diagnostics_file, tracers, group)

        end

    else

        dimensions!(diagnostics_file, computed_output)
        save_computed_output!(diagnostics_file, computed_output, group)
        potential_and_background_potential_energy!(diagnostics_file, computed_output, tracers, group)
        save_snaphots!(diagnostics_file, tracers, velocities, group; snapshots = 1:25)
        φ_interface_flux!(diagnostics_file, tracers, :S, group)
        φ_interface_flux!(diagnostics_file, tracers, :T, group)
        ha_φ_flux!(diagnostics_file, tracers, :S_ha, group)
        ha_φ_flux!(diagnostics_file, tracers, :T_ha, group)
        ∫gρw!(diagnostics_file, computed_output, velocities, group)
        initial_non_dim_numbers!(diagnostics_file, computed_output, group)
        save_horizontally_averaged_fields!(diagnostics_file, computed_output, tracers, group)
        compute_Ẽ!(diagnostics_file, computed_output, tracers, group, interface_offset)
        interface_thickness!(diagnostics_file, tracers, group)

    end

    return nothing
end
"""
    function save_diagnostic!(diagnostics_file::AbstractString, diagnostic_function!::Function,
                                function_args)
Save `diagnostic_function!` to `diagnostics_file`. The `function_args` are the arguments for
`diagnostic_function`. This allows a nice way to add other diagnostics that follows syntax
I have used above though this is really a nice to have rather than necessary function.
`function_args` will always require the data to compute the diagnostic from, the `group` to
save it to in `diagnostics_file`. **Note:** the `group` must end with `\` i.e. `lineareos\`
will create the group in `diagnostics_file` but `lineareos` will just prepend variable name.
Either will work but best to save in groups.
"""
function save_diagnostic!(diagnostics_file::AbstractString, diagnostic_function!::Function,
                            function_args)

    diagnostic_function!(diagnostics_file, function_args...)
    return nothing
end
"Local Kolmogorov length scale."
η(ν, ε) = (ν^3 / ε)^(1/4)
"Local Batchelor length scale."
Ba(η, Sc) = η / sqrt(Sc)
"""
    function save_computed_output!(diagnostics_file::AbstractString, computed_output::AbstractString, group)
Save diagnostics from `computed_output` (energetics, Rᵨ, horizontally averaged N² and density)
"""
function save_computed_output!(diagnostics_file::AbstractString, computed_output::AbstractString, group)

    NCDataset(computed_output) do ds

        if isfile(diagnostics_file)
            jldopen(diagnostics_file, "a+") do file
                file[group*"R_ρ"] = ds["R_ρ"][:]
                file[group*"∫ε"] = ds["∫ε"][:]
                file[group*"∫Eₖ"] = ds["∫Eₖ"][:]
                file[group*"∫wb"] = ds["∫wb"][:]
                ν, Sc = ds.attrib["ν (m²s⁻¹)"], ds.attrib["Sc"]
                η_ = η.(ν, ds["ε_maximum"][:])
                file[group*"η"] = η_
                file[group*"Ba"] = Ba.(η_, Sc)
            end
        else
            jldopen(diagnostics_file, "w") do file
                file[group*"R_ρ"] = ds["R_ρ"][:]
                file[group*"∫ε"] = ds["∫ε"][:]
                file[group*"∫Eₖ"] = ds["∫Eₖ"][:]
                file[group*"∫wb"] = ds["∫wb"][:]
                ν, Sc = ds.attrib["ν (m²s⁻¹)"], ds.attrib["Sc"]
                η_ = η.(ν, ds["∫ε"][:])
                file[group*"η"] = η_
                file[group*"Ba"] = Ba(η_, Sc)
            end
        end

    end

    return nothing
end
"""
    function dimensions!(diagnostics_file::AbstractString, co::AbstractString)
Save space, time and any other variable that acts as a dimension (e.g `z✶`).
"""
function dimensions!(diagnostics_file::AbstractString, co::AbstractString)

    dims = ("time", "x_caa", "y_aca", "z_aac", "z_aaf")
    NCDataset(co) do ds

        if isfile(diagnostics_file)
            jldopen(diagnostics_file, "a+") do file
                for d ∈ dims
                    file["dims/"*d] = ds[d][:]
                end
            end
        else
            jldopen(diagnostics_file, "w") do file
                for d ∈ dims
                    file["dims/"*d] = ds[d][:]
                end
            end
        end

    end

    return nothing
end

function initial_non_dim_numbers!(diagnostics_file::AbstractString, computed_output::AbstractString, group)

    nd_nums = ("Pr", "Sc", "Le", "RaS", "RaT", "Rᵨ", )
    other_attribs = ("interface_depth", "EOS", "Reference density (kgm⁻³)", "ν (m²s⁻¹)",
                     "κₛ (m²s⁻¹)", "κₜ (m²s⁻¹)")

    NCDataset(computed_output) do ds

        if isfile(diagnostics_file)
            jldopen(diagnostics_file, "a+") do file
                for nd ∈ nd_nums
                    file[group*"attrib/"*nd] = ds.attrib[nd]
                end
                for oa ∈ other_attribs
                    file[group*"attrib/"*oa] = ds.attrib[oa]
                end
            end
        else
            jldopen(diagnostics_file, "w") do file
                for nd ∈ nd_nums
                    file[group*"attrib/"*nd] = ds.attrib[nd]
                end
                for oa ∈ other_attribs
                    file[group*"attrib/"*oa] = ds.attrib[oa]
                end
            end
        end

    end
    return nothing
end
"""
    function update_diagnostic!(diagnostics_file::AbstractString, group::AbstractString,
                                key::AbstractString, tracers::AbstractString,
                                computed_output::AbstractString)
Update the diagnostic (i.e. run again) at `key` in diagnostics_file. This function first deletes
what is saved at `key` then calls [save_diagnostics](@ref) which should only update the relevant
diagnostic. **Note:** the function will find the relevant `key` then remove all keys that the
function saves so they can all be recomputed.
"""
function update_diagnostic!(diagnostics_file::AbstractString, group::AbstractString,
                            key::AbstractString, tracers::AbstractString,
                            computed_output::AbstractString;
                            interface_offset = 4)

    S_flux_keys = ("S_flux", "S_interface_idx")
    ha_S_flux_keys = ("ha_S_flux", "ha_S_interface_idx")
    T_flux_keys = ("T_flux", "T_interface_idx")
    ha_T_flux_keys = ("ha_T_flux", "ha_T_interface_idx")
    interface_thickness_keys = ("hₜ", "hₛ", "r","ΔS", "ΔT")
    Ẽ_keys = ("Ẽ", "Tₗ_Tᵤ_ts", "Sₗ_Sᵤ_ts", "ρₗ_ρᵤ_ts")

      keys_to_remove =  if key ∈ S_flux_keys
                            S_flux_keys
                        elseif key ∈ T_flux_keys
                            T_flux_keys
                        elseif key ∈ ha_S_flux_keys
                            ha_S_flux_keys
                        elseif key ∈ ha_T_flux_keys
                            ha_T_flux_keys
                        elseif key ∈ interface_thickness_keys
                            interface_thickness_keys
                        elseif key ∈ Ẽ_keys
                            Ẽ_keys
                        end

    delete_keys!(diagnostics_file, group, keys_to_remove)

    group = group * '/'
    if keys_to_remove == S_flux_keys
        φ_interface_flux!(diagnostics_file, tracers, :S, group)
    elseif keys_to_remove == T_flux_keys
        φ_interface_flux!(diagnostics_file, tracers, :T, group)
    elseif keys_to_remove == ha_S_flux_keys
        ha_φ_flux!(diagnostics_file, tracers, :S, group)
    elseif keys_to_remove == ha_T_flux_keys
        ha_φ_flux!(diagnostics_file, tracers, :T, group)
    elseif keys_to_remove == Ẽ_keys
        compute_Ẽ!(diagnostics_file, computed_output, tracers, group, interface_offset)
    elseif keys_to_remove == interface_thickness_keys
        interface_thickness!(diagnostics_file, tracers, group)
    end

    return nothing
end
function delete_keys!(diagnostics_file, group, keys_to_remove)

    jldopen(diagnostics_file, "a+") do f
        for k ∈ keys_to_remove
            delete!(f, group*"/"*k)
        end
    end

    return nothing
end
"""
    function φ_interface_flux!(diagnostics_file::AbstractString, tracers::AbstractString, tracer::Symbol)
Calculate the flux through the diffusive interface for the `φ` tracer. The diffusive interface
is found in the sorted profile where φ > Δφ / 2 (with the previous two levels saved so it is
not a single value). The flux is then calculated as the change in φ content to
z✶(Δφ / 2) were z✶ is the backgrund z profile.
"""
function φ_interface_flux!(diagnostics_file::AbstractString, tracers::AbstractString, tracer::Symbol, group)

    NCDataset(tracers) do ds

        φ = ds[tracer]
        Δφ₀ = φ[1, 1, 1, 1] - 0.5 * (φ[1, 1, 1, 1] - φ[1, 1, end, 1])
        timestamps = ds[:time][:]
        Δt = diff(timestamps)
        Δx, Δy, Δz = ds[:Δx_caa][1], ds[:Δy_aca][1], ds[:Δz_aac][1]
        ΔV = Δx * Δy * Δz
        V = (1:length(reshape(φ[:, :, :, 1], :))) * ΔV
        Lx, Ly = Δx * ds.dim[:x_caa], Δy * ds.dim[:y_aca]
        SA = Lx * Ly
        z✶ = V / SA
        save_z✶!(diagnostics_file, z✶)
        Δz✶ = diff(z✶)[1]

        φ_interface_flux = Array{Float64}(undef, 3, length(Δt))
        interface_idx = Array{Int64}(undef, length(Δt))

        for i ∈ eachindex(Δt)

            φₜ = [reshape(φ[:, :, :, i], :) reshape(φ[:, :, :, i+1], :)]
            sort!(φₜ, rev = true, dims = 1)
            ∫φdz✶ = cumsum(φₜ * Δz✶, dims = 1)
            dₜ∫φdz✶ = vec(diff(∫φdz✶, dims = 2) ./ Δt[i])

            interface_idx[i] = ii = findfirst(φₜ[:, 1].< Δφ₀) - 1
            interface_idxs = [ii-1, ii, ii+1]

            φ_interface_flux[:, i] .= dₜ∫φdz✶[interface_idxs]

        end

        save_fluxes!(diagnostics_file, φ_interface_flux, interface_idx, tracer, group)

    end

    return nothing
end
"""
    function ha_φ_flux!(diagnostics_file::AbstractString, tracers::AbstractString, tracer::Symbol)
Calculate the flux through the diffusive interface for the horizontally averaged `φ` tracer.
The diffusive interface is found in the horizontally averaged, sorted profile where φ > Δφ / 2
(with the previous two levels saved so it is not a single value). The flux is then calculated as the change in φ content to
z(Δφ / 2) were z is the depth. In this case the whole flux field is saved and the index of
the interface is returned separately.
"""
function ha_φ_flux!(diagnostics_file::AbstractString, tracers::AbstractString, tracer::Symbol, group)

    NCDataset(tracers) do ds

        φ = ds[tracer]
        Δφ₀ = φ[1, 1] - 0.5 * (φ[1, 1] - φ[end, 1])
        timestamps = ds[:time][:]
        Δt = diff(timestamps)
        z = ds[:z_aac][:]
        Δz = ds[:Δz_aac][1]

        φ_flux = Array{Float64}(undef, length(z), length(Δt))
        interface_idx = Array{Int64}(undef, length(Δt))

        for i ∈ eachindex(Δt)

            φₜ = [φ[:, i] φ[:, i+1]]
            sort!(φₜ, rev = true, dims = 1)
            ∫φdz = cumsum(φₜ * Δz, dims = 1)
            dₜ∫φdz = vec(diff(∫φdz, dims = 2) ./ Δt[i])

            interface_idx[i] = findfirst(φₜ[:, 1] .< Δφ₀) - 1

            φ_flux[:, i] .= dₜ∫φdz

        end

        ha_tracer = string(tracer)
        save_fluxes!(diagnostics_file, φ_flux, interface_idx, ha_tracer, group)

    end

    return nothing
end
"""
    function save_fluxes!(diagnostics_file, φ_interface_flux, interface_idx, tracer)
Save the flux through, and index of, the diffusive interface. **Note** the index is for the reshaped
and resorted vector.
"""
function save_fluxes!(diagnostics_file, φ_interface_flux, interface_idx, tracer, group)

    if isfile(diagnostics_file)
        jldopen(diagnostics_file, "a+") do file
            file[group*string(tracer)*"_flux"] = φ_interface_flux
            file[group*string(tracer)*"_interface_idx"] = interface_idx
        end
    else
        jldopen(diagnostics_file, "w") do file
            file[group*string(tracer)*"_flux"] = φ_interface_flux
            file[group*string(tracer)*"_interface_idx"] = interface_idx
        end
    end

    return nothing
end
function save_z✶!(diagnostics_file, z✶)

    f = jldopen(diagnostics_file, "a+")
    if !haskey(f, "dims/z✶")
        f["dims/z✶"] =  z✶
    end
    close(f)
    return nothing
end
"""
    function interface_thickness!(diagnostics_file::AbstractString, tracers::AbstractString)
Calculate the interface thickness using the method in appendix A of
[Sommer et al. (2013)](https://journals.ametsoc.org/view/journals/atot/30/8/jtech-d-12-00272_1.xml).
"""
function interface_thickness!(diagnostics_file::AbstractString, tracers::AbstractString, group)

    ds = NCDataset(tracers) do ds

        timestamps = ds[:time][:]

        zC = abs.(reverse(ds[:z_aac][:]))

        Tmidpoint = 0.5 * (ds[:T_ha][1] .+ ds[:T_ha][end])
        Smidpoint = 0.5 * (ds[:S_ha][1] .+ ds[:S_ha][end])

        hₜ = similar(timestamps[1:end-1])
        hₛ = similar(timestamps[1:end-1])
        ΔT = similar(timestamps[1:end-1])
        ΔS = similar(timestamps[1:end-1])

        interface_offfset = 50 # how far from the interface the mean S and T for a layer are calculate

        for t ∈ eachindex(hₜ)
            T = ds[:T_ha][:, t+1]
            find_interface = findfirst(T .< Tmidpoint)
            ΔT[t] = abs.(mean(ds[:T_ha][1:find_interface-interface_offfset]) .- mean(ds[:T_ha][find_interface+interface_offfset:end]))
            find = findall(Tmidpoint - (ΔT[t]/4) .<  T .< Tmidpoint + (ΔT[t]/4))
            intercept, slope = [ones(length(find)) zC[find]] \ T[find]
            hₜ[t] = ΔT[t] / slope
        end

        for t ∈ eachindex(hₛ)
            S = ds[:S_ha][:, t+1]
            find_interface = findfirst(S .< Smidpoint)
            ΔS[t] = abs.(mean(ds[:S_ha][1:find_interface-interface_offfset]) .- mean(ds[:S_ha][find_interface+interface_offfset:end]))
            find = findall(Smidpoint - (ΔS[t]/4) .<  S .< Smidpoint + (ΔS[t]/4))
            intercept, slope = [ones(length(find)) zC[find]] \ S[find]
            hₛ[t] = ΔS[t] / slope
        end

        r = hₜ ./ hₛ
        save_interface_thickness!(diagnostics_file, hₜ, hₛ, r, ΔS, ΔT, group)
    end

    return nothing
end
function save_interface_thickness!(diagnostics_file, hₜ, hₛ, r, ΔS, ΔT, group)

    if isfile(diagnostics_file)
        jldopen(diagnostics_file, "a+") do file
            file[group*"hₜ"] = hₜ
            file[group*"hₛ"] = hₛ
            file[group*"r"] = r
            file[group*"ΔS"] = ΔS
            file[group*"ΔT"] = ΔT
        end
    else
        jldopen(diagnostics_file, "w") do file
            file[group*"hₜ"] = hₜ
            file[group*"hₛ"] = hₛ
            file[group*"r"] = r
            file[group*"ΔS"] = ΔS
            file[group*"ΔT"] = ΔT
        end
    end

    return nothing
end
"""
    function compute_Ẽ!(diagnostics_file::AbstractString, co::AbstractString, tracers::AbstractString, group;
               interface_offset = 1)
Compute the entrainment parameter using equation (y) from [McDougall (1981)](https://www.sciencedirect.com/science/article/pii/007966118190001X).
This is defined as
```math
Ẽ = \\frac{ρₗHₗdTₗ/dt + ρᵤHᵤdTᵤ/dt}{ρₗHₗdTₗ/dt - ρᵤHᵤdTᵤ/dt}.
```
The interface is found as the midpoint temperature and salinity of the two layers and the
heights Hᵤ and Hₗ are calculated from this but they take initial values of the domain divided
in two as that is how the experiments are set.
Also saved is a time series of the average temperature, salinity and density in the upper and
lower layers.
The warg `interface_offset` is how far away to move (vertically) from the interface to define
the upper and lower layers. This is passed by [save_diagnostics!](@ref)
"""
function compute_Ẽ!(diagnostics_file::AbstractString, co::AbstractString, tracers::AbstractString,
                    group, interface_offset)

    ds_co = NCDataset(co)
    ds_tracers = NCDataset(tracers)

    timestamps = ds_tracers[:time][:]
    zC = ds_tracers[:z_aac][:]
    Δx, Δy, Δz = ds_tracers[:Δx_caa][1], ds_tracers[:Δy_aca][1], ds_tracers[:Δz_aac][1]

    Δt = diff(timestamps)
    Ẽ = similar(Δt)
    T̄_ts = Array{Float64}(undef, length(Δt), 2)
    S̄_ts = similar(T̄_ts)
    ρ̄_ts = similar(T̄_ts)

    S = ds_tracers[:S_ha]
    T = ds_tracers[:T_ha]
    mid_T = 0.5 * (T[end, 1] + T[1, 1])
    σ = ds_co[:σ_ha]

    for t ∈ eachindex(Δt)

        T_ha = T[:, t]
        interface_idx = findfirst(T_ha .< mid_T)
        lower = 1:interface_idx-1-interface_offset # the extra -1 is needed because of `findfirst` index.
        upper = interface_idx+interface_offset:length(zC)

        Hₗ = abs(zC[lower][1] - zC[lower][end])
        Hᵤ = abs(zC[upper][1] - zC[upper][end])

        Tₗ = T[lower, t:t+1]
        T̄ₗ = mean(Tₗ, dims = 1)
        dₜT̄ₗ = diff(T̄ₗ, dims = 2) ./ Δt[t]
        _ρₗ = 0.5 * sum(σ[lower, t:t+1], dims = 2)
        ρₗ = mean(_ρₗ)

        Tᵤ = T[upper, t:t+1]
        T̄ᵤ = mean(Tᵤ, dims = 1)
        dₜT̄ᵤ = diff(T̄ᵤ, dims = 2) ./ Δt[t]
        _ρᵤ = 0.5 * sum(σ[upper, t:t+1], dims = 2)
        ρᵤ = mean(_ρᵤ)

        S̄ₗ, S̄ᵤ = mean(S[lower, t:t+1], dims = 1), mean(S[upper, t:t+1], dims = 1)
        S̄_ts[t, :] .= [S̄ₗ[1, 1], S̄ᵤ[1, 1]]
        T̄_ts[t, :] .= [T̄ₗ[1, 1], T̄ᵤ[1, 1]]
        ρ̄_ts[t, :] .= [ρₗ, ρᵤ]

        _Ẽ = @. (ρₗ * Hₗ * dₜT̄ₗ + ρᵤ * Hᵤ * dₜT̄ᵤ) / (ρₗ * Hₗ * dₜT̄ₗ - ρᵤ * Hᵤ * dₜT̄ᵤ)
        Ẽ[t] = _Ẽ[1, 1]
    end

    close(ds_co)
    close(ds_tracers)

    save_Ẽ!(diagnostics_file, Ẽ, T̄_ts, S̄_ts, ρ̄_ts, group)

    return nothing
end
function save_Ẽ!(diagnostics_file, Ẽ, T̄_ts, S̄_ts, ρ̄_ts, group)

    if isfile(diagnostics_file)
        jldopen(diagnostics_file, "a+") do file
            file[group*"Ẽ"] = Ẽ
            file[group*"Tₗ_Tᵤ_ts"] = T̄_ts
            file[group*"Sₗ_Sᵤ_ts"] = S̄_ts
            file[group*"ρₗ_ρᵤ_ts"] = ρ̄_ts
        end
    else
        jldopen(diagnostics_file, "w") do file
            file[group*"Ẽ"] = Ẽ
            file[group*"Tₗ_Tᵤ_ts"] = T̄_ts
            file[group*"Sₗ_Sᵤ_ts"] = S̄_ts
            file[group*"ρₗ_ρᵤ_ts"] = ρ̄_ts
        end
    end
    return nothing
end
"""
    function ∫gρw!(diagnostics_file::AbstractString, computed_output::AbstractString, velocities::AbstractString)
Compute the buoyancy flux from model density and vertical velocity fields.
"""
function ∫gρw!(diagnostics_file::AbstractString, computed_output::AbstractString, velocities::AbstractString, group)

    ds_co = NCDataset(computed_output, "a")
    timestamps = ds_co[:time][:]
    ρ₀ = ds_co.attrib["Reference density (kgm⁻³)"]
    ΔV = ds_co[:Δx_caa][1] * ds_co[:Δy_aca][1] * ds_co[:Δz_aac][1]
    ds_vel = NCDataset(velocities)

    g = 9.81
    ∫gρw = similar(timestamps)
    for t ∈ eachindex(timestamps)

        σ = ds_co[:σ][:, :, :, t]
        σ1 = @view σ[:, :, 1:end-1]
        σ2 = @view σ[:, :, 2:end]
        σ_interp = cat(σ[:, :, 1], 0.5 * (σ1 .+ σ2), σ[:, :, end], dims = 3)
        w = ds_vel[:w][:, :, :, t]

        ∫gρw[t] = (g / ρ₀) * sum(σ_interp .* w) * ΔV

    end
    close(ds_co)
    close(ds_vel)

    save_∫gρw!(diagnostics_file, ∫gρw, group)

    return nothing
end
function save_∫gρw!(diagnostics_file, ∫gρw, group)

    if isfile(diagnostics_file)
        jldopen(diagnostics_file, "a+") do file
            file[group*"∫gρw"] = ∫gρw
        end
    else
        jldopen(diagnostics_file, "w") do file
            file[group*"∫gρw"] = ∫gρw
        end
    end

    return nothing
end
function save_horizontally_averaged_fields!(diagnostics_file::AbstractString,
                                            computed_output::AbstractString,
                                            tracers::AbstractString, group)
    NCDataset(tracers) do ds
        if isfile(diagnostics_file)
            jldopen(diagnostics_file, "a+") do file
                file[group*"T_ha"] = ds[:T_ha][:, :]
                file[group*"S_ha"] = ds[:S_ha][:, :]
            end
        else
            jldopen(diagnostics_file, "w") do file
                file[group*"T_ha"] = ds[:T_ha][:, :]
                file[group*"S_ha"] = ds[:S_ha][:, :]
            end
        end
    end

    NCDataset(computed_output) do ds
        if isfile(diagnostics_file)
            jldopen(diagnostics_file, "a+") do file
                file[group*"σ_ha"] = ds[:σ_ha][:, :]
                file[group*"N²_ha"] = ds[:N²_ha][:, :]
            end
        else
            jldopen(diagnostics_file, "w") do file
                file[group*"σ_ha"] = ds[:σ_ha][:, :]
                file[group*"N²_ha"] = ds[:N²_ha][:, :]
            end
        end
    end

    return nothing
end
"""
    function potential_and_background_potential_energy!(computed_output::AbstractString)
Compute and append the potential and background energy to `computed_output`. **Note** the
PE and BPE are both referenced to ``z = 0``, and the saved quantities are volume integrated.
"""
function potential_and_background_potential_energy!(diagnostics_file::AbstractString,
                                                    computed_output::AbstractString,
                                                    tracers::AbstractString, group)

    NCDataset(computed_output, "a") do ds

        SA = sum(ds["Δx_caa"][:]) * sum(ds["Δy_aca"][:])
        t = ds[:time][:]
        ρ₀ = ds.attrib["Reference density (kgm⁻³)"]
        ΔV = ds[:Δx_caa][1] * ds[:Δy_aca][1] * ds[:Δz_aac][1]
        σ = ds[:σ]

        V = cumsum(ones(length(reshape(σ[:, :, :, 1], :)))) * ΔV
        z✶ = V / SA

        ds_tracers = NCDataset(tracers)
        T = ds_tracers[:T]
        mid_T = 0.5 * (T[1, 1, end, 1] + T[1, 1, 1, 1])

        # background potential energy
        Eb = similar(t)
        Eb_lower = similar(t)
        Eb_upper = similar(t)
        T_interface_z✶ = similar(t)
        g = 9.81
        for i ∈ eachindex(t)
            σᵢ = σ[:, :, :, i] .- ρ₀
            σᵢ_array = reshape(σᵢ, :)
            sort!(σᵢ_array, rev = true)
            Eb[i] = (g / ρ₀) * sum(σᵢ_array .* z✶ * ΔV)
            # find index of interface
            Tᵢ = reshape(T[:, :, :, i], :)
            sort!(Tᵢ, rev = true)
            T_interface = findfirst(Tᵢ .≤ mid_T) - 1
            T_interface_z✶[i] = z✶[T_interface]
            # compute BPE in each layer
            Eb_lower[i] = (g / ρ₀) * sum(σᵢ_array[1:T_interface] .* z✶[1:T_interface] * ΔV)
            Eb_upper[i] = (g / ρ₀) * sum(σᵢ_array[T_interface+1:end] .* z✶[T_interface+1:end] * ΔV)
        end

        # potential energy
        z = ds["z_aac"]
        Nx = ds.group["grid_reconstruction"].attrib["Nx"]
        Ny = ds.group["grid_reconstruction"].attrib["Ny"]
        Nz = ds.group["grid_reconstruction"].attrib["Nz"]
        z_ref0 = reverse(abs.(z))
        z_grid = reshape(repeat(z_ref0, inner = Nx * Ny), (Nx, Ny, Nz))
        Ep = similar(t)
        Ep_lower = similar(t)
        Ep_upper = similar(t)
        for i ∈ eachindex(t)
            σᵢ = σ[:, :, :, i] .- ρ₀
            Ep[i] = (g / ρ₀) * sum(σᵢ .* z_grid * ΔV)
            # find index of interface
            T_interface = findfirst(z_ref0 .≥ T_interface_z✶[i]) - 1
            # compute PE within each layer
            Ep_lower[i] = (g / ρ₀) * sum(σᵢ[:, :, 1:T_interface] .* z_grid[:, :, 1:T_interface] * ΔV)
            Ep_upper[i] = (g / ρ₀) * sum(σᵢ[:, :, T_interface+1:end] .* z_grid[:, :, T_interface+1:end] * ΔV)
        end

    close(ds_tracers)

    save_pe_and_bpe!(diagnostics_file, Ep, Ep_lower, Ep_upper, Eb, Eb_lower, Eb_upper, group)

    end

    return nothing
end
function save_pe_and_bpe!(diagnostics_file::AbstractString, Ep, Ep_lower, Ep_upper, Eb, Eb_lower, Eb_upper, group)

    if isfile(diagnostics_file)
        jldopen(diagnostics_file, "a+") do file
            file[group*"Ep"] = Ep
            file[group*"Ep_lower"] = Ep_lower
            file[group*"Ep_upper"] = Ep_upper
            file[group*"Eb"] = Eb
            file[group*"Eb_lower"] = Eb_lower
            file[group*"Eb_upper"] = Eb_upper
        end
    else
        jldopen(diagnostics_file, "w") do file
            file[group*"Ep"] = Ep
            file[group*"Ep_lower"] = Ep_lower
            file[group*"Ep_upper"] = Ep_upper
            file[group*"Eb"] = Eb
            file[group*"Eb_lower"] = Eb_lower
            file[group*"Eb_upper"] = Eb_upper
        end
    end
    return nothing
end
function save_snaphots!(diagnostics_file::AbstractString, tracers::AbstractString,
                        velocities::AbstractString, group; snapshots = 1:25)

    if isfile(diagnostics_file)

        jldopen(diagnostics_file, "w") do file
            NCDataset(tracers) do ds

                t = ds["time"][:]
                for i ∈ snapshots
                    file[group*"S/xzslice_$(t[i])"] = ds[:S][:, 1, :, i]
                    file[group*"T/yzslice_$(t[i])"] = ds[:T][1, :, :, i]
                end
            end

            NCDataset(velocities) do ds

                t = ds["time"][:]
                for i ∈ snapshots
                    file[group*"w/w_zmean_$(t[i])"] = mean(ds[:w][:, :, :, i], dims = 3)
                end
            end
        end

    else

        jldopen(diagnostics_file, "w") do file
            NCDataset(tracers) do ds

                t = ds["time"][:]
                for i ∈ snapshots
                    file["S/xzslice_$(t[i])"] = ds[:S][:, 1, :, i]
                    file["T/yzslice_$(t[i])"] = ds[:T][1, :, :, i]
                end
            end

            NCDataset(velocities) do ds

                t = ds["time"][:]
                for i ∈ snapshots
                    file["w/w_zmean_$(t[i])"] = mean(ds[:w][:, :, :, i], dims = 3)
                end
            end
        end

    end

    return nothing
end
