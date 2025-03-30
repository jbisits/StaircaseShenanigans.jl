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

    interface_depth = NCDataset(computed_output) do co
                          co.attrib[:interface_depth]
                      end
    ds = NCDataset(tracers)

    S_u = S_g = ds[:Sᵤ_mean][:]
    S_l = S_f = ds[:Sₗ_mean][:]
    T_u = T_f = ds[:Tᵤ_mean][:]
    T_l = T_g = ds[:Tₗ_mean][:]

    eos_vec = fill(eos, length(S_u))
    ρ_u = total_density.(T_u, S_u, interface_depth, eos_vec)
    ρ_l = total_density.(T_l, S_l, interface_depth, eos_vec)
    ρ_f = total_density.(T_f, S_f, interface_depth, eos_vec)
    ρ_g = total_density.(T_g, S_g, interface_depth, eos_vec)

    R_ρ = @. (0.5 * (ρ_f - ρ_u) + 0.5 * (ρ_l - ρ_g)) / (0.5 * (ρ_f - ρ_l) + 0.5 * (ρ_u - ρ_g))

    NCDataset(computed_output, "a") do ds2
        if haskey(ds2, "R_ρ")
            # rename so if picking up can saved the whole array
            renameVar(ds2, "R_ρ", "R_ρ_prior_cp")
            defVar(ds2, "R_ρ", R_ρ, ("time",),
                    attrib = Dict("long_name" => "Density ratio"))
        else
            defVar(ds2, "R_ρ", R_ρ, ("time",),
                    attrib = Dict("long_name" => "Density ratio"))
        end
    end

    close(ds)

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
                           computed_output::AbstractString;
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
            φ_interface_flux!(diagnostics_file, tracers, :S, group)
            φ_interface_flux!(diagnostics_file, tracers, :T, group)
            ha_φ_flux!(diagnostics_file, tracers, :S, group)
            ha_φ_flux!(diagnostics_file, tracers, :T, group)
            compute_Ẽ!(diagnostics_file, computed_output, tracers, group, interface_offset)
            interface_thickness!(diagnostics_file, tracers, group)
        end

    else

        dimensions!(diagnostics_file, computed_output)
        save_computed_output!(diagnostics_file, computed_output, group)
        φ_interface_flux!(diagnostics_file, tracers, :S, group)
        φ_interface_flux!(diagnostics_file, tracers, :T, group)
        ha_φ_flux!(diagnostics_file, tracers, :S, group)
        ha_φ_flux!(diagnostics_file, tracers, :T, group)
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

    dims = ("time", "x_caa", "y_aca", "z_aaf")
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

            φₜ_interp = 0.5 * vec(sum(φₜ, dims = 2))
            interface_idx[i] = ii = findfirst(φₜ_interp .< Δφ₀) - 1
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
        Δφ₀ = φ[1, 1, 1, 1] - 0.5 * (φ[1, 1, 1, 1] - φ[1, 1, end, 1])
        timestamps = ds[:time][:]
        Δt = diff(timestamps)
        z = ds[:z_aac][:]
        Δz = ds[:Δz_aac][1]

        φ_flux = Array{Float64}(undef, length(z), length(Δt))
        interface_idx = Array{Int64}(undef, length(Δt))

        for i ∈ eachindex(Δt)

            φₜ = [reshape(mean(φ[:, :, :, i], dims = (1, 2)), :) reshape(mean(φ[:, :, :, i+1], dims = (1, 2)), :)]
            sort!(φₜ, rev = true, dims = 1)
            ∫φdz = cumsum(φₜ * Δz, dims = 1)
            dₜ∫φdz = vec(diff(∫φdz, dims = 2) ./ Δt[i])

            φₜ_interp = 0.5 * vec(sum(φₜ, dims = 2))
            interface_idx[i] = findfirst(φₜ_interp .< Δφ₀) - 1

            φ_flux[:, i] .= dₜ∫φdz

        end

        ha_tracer = "ha_"*string(tracer)
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
        ΔT = abs.(ds[:Tₗ_mean][:] .- ds[:Tᵤ_mean][:])
        ΔS = abs.(ds[:Sₗ_mean][:] .- ds[:Sᵤ_mean][:])
        zC = abs.(ds[:z_aac][:])

        Tmidpoint = 0.5 * (ds[:Tₗ_mean][1] .+ ds[:Tᵤ_mean][1])
        Smidpoint = 0.5 * (ds[:Sₗ_mean][1] .+ ds[:Sᵤ_mean][1])

        hₜ = similar(timestamps[1:end-1])
        hₛ = similar(timestamps[1:end-1])

        for t ∈ eachindex(hₜ)
            T = reshape(mean(ds[:T][:, :, :, t+1], dims = (1, 2)), :)
            find = findall(Tmidpoint - (ΔT[t+1]/4) .<  T .< Tmidpoint + (ΔT[t+1]/4))
            intercept, slope = [ones(length(find)) zC[find]] \ T[find]
            hₜ[t] = ΔT[t] / slope
        end

        for t ∈ eachindex(hₛ)
            S = reshape(mean(ds[:S][:, :, :, t+1], dims = (1, 2)), :)
            find = findall(Smidpoint - (ΔS[t+1]/4) .<  S .< Smidpoint + (ΔS[t+1]/4))
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

    timestamps = ds_co[:time][:]
    zC = ds_co[:z_aac][:]
    Δx, Δy, Δz = ds_co[:Δx_caa][1], ds_co[:Δy_aca][1], ds_co[:Δz_aac][1]

    Δt = diff(timestamps)
    Ẽ = similar(Δt)
    T̄_ts = Array{Float64}(undef, length(Δt), 2)
    S̄_ts = similar(T̄_ts)
    ρ̄_ts = similar(T̄_ts)

    S = ds_tracers[:S]
    T = ds_tracers[:T]
    mid_T = 0.5 * (T[1, 1, end, 1] + T[1, 1, 1, 1])
    σ = ds_co[:σ]

    for t ∈ eachindex(Δt)

        T_ha = reshape(mean(T[:, :, :, t], dims = (1, 2)), :)
        interface_idx = findfirst(T_ha .< mid_T)
        lower = 1:interface_idx-1-interface_offset # the extra -1 is needed because of `findfirst` index.
        upper = interface_idx+interface_offset:length(zC)

        Hₗ = abs(zC[lower][1] - zC[lower][end])
        Hᵤ = abs(zC[upper][1] - zC[upper][end])

        Tₗ = T[:, :, lower, t:t+1]
        T̄ₗ = mean(Tₗ, dims = (1, 2, 3))
        dₜT̄ₗ = diff(T̄ₗ, dims = 4) ./ Δt[t]
        _ρₗ = 0.5 * sum(σ[:, :, lower, t:t+1], dims = 4)
        ρₗ = mean(_ρₗ)

        Tᵤ = T[:, :, upper, t:t+1]
        T̄ᵤ = mean(Tᵤ, dims = (1, 2, 3))
        dₜT̄ᵤ = diff(T̄ᵤ, dims = 4) ./ Δt[t]
        _ρᵤ = 0.5 * sum(σ[:, :, upper, t:t+1], dims = 4)
        ρᵤ = mean(_ρᵤ)

        S̄ₗ, S̄ᵤ = mean(S[:, :, lower, t:t+1], dims = (1, 2, 3)), mean(S[:, :, upper, t:t+1], dims = (1, 2, 3))
        S̄_ts[t, :] .= [S̄ₗ[1, 1, 1, 1], S̄ᵤ[1, 1, 1, 1]]
        T̄_ts[t, :] .= [T̄ₗ[1, 1, 1, 1], T̄ᵤ[1, 1, 1, 1]]
        ρ̄_ts[t, :] .= [ρₗ, ρᵤ]

        _Ẽ = @. (ρₗ * Hₗ * dₜT̄ₗ + ρᵤ * Hᵤ * dₜT̄ᵤ) / (ρₗ * Hₗ * dₜT̄ₗ - ρᵤ * Hᵤ * dₜT̄ᵤ)
        Ẽ[t] = _Ẽ[1, 1, 1, 1]
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
