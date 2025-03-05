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
    function save_diagnostics!(diagnostics_file::AbstractString, tracers::AbstractString, computed_output::AbstractString)
Save the diagnostics in this script, diagnostics.jl, to `diagnostics_file`.
"""
function save_diagnostics!(diagnostics_file::AbstractString, tracers::AbstractString,
                           computed_output::AbstractString; eos = nothing)

    eos = isnothing(eos) ? "" : eos[end] == '/' ? eos : eos * "/" # creates a group in the saved output.
    if isfile(diagnostics_file)

        haskey(diagnostics_file, eos*"S_flux") ? nothing : φ_interface_flux!(diagnostics_file, tracers, :S, eos)
        haskey(diagnostics_file, eos*"T_flux") ? nothing : φ_interface_flux!(diagnostics_file, tracers, :T, eos)
             haskey(diagnostics_file, eos*"Ẽ") ? nothing : Ẽ!(diagnostics_file, computed_output, eos)
            haskey(diagnostics_file, eos*"hₜ") ? nothing : interface_thickness!(diagnostics_file, tracers, eos)

    else

        φ_interface_flux!(diagnostics_file, tracers, :S, eos)
        φ_interface_flux!(diagnostics_file, tracers, :T, eos)
        Ẽ!(diagnostics_file, computed_output, eos)
        interface_thickness!(diagnostics_file, tracers, eos)

    end

    return nothing
end

"""
    function update_diagnostic!(diagnostics_file::AbstractString, key::AbstractString, tracers::AbstractString,
                                computed_output::AbstractString; eos = nothing)
Update the diagnostic (i.e. run again) at `key` in diagnostics_file. This function first deletes
what is saved at `key` then calls [save_diagnostics](@ref) which should only update the relevant
diagnostic. **Note:** the function will find the relevant `key` then remove all keys that the
function saves so they can all be recomputed
"""
function update_diagnostic!(diagnostics_file::AbstractString, key::AbstractString, tracers::AbstractString,
                            computed_output::AbstractString; eos = nothing)

    eos = isnothing(eos) : "" : eos * "/"
    S_flux_keys = eos .* string.([S_flux, S_interface_idx])
    interface_thickness_keys = eos .* string.([hₜ, hₛ, r, ΔS, ΔT])
    Ẽ_keys = eos .* string.([Ẽ, Tₗ_Tᵤ_ts, Sₗ_Sᵤ_ts, ρₗ_ρᵤ_ts])

    keys =  if key ∈ S_flux_keys
                S_flux_keys
            elseif key ∈ interface_thickness_keys
                interface_thickness_keys
            elseif key ∈ Ẽ_keys
                Ẽ_keys
            end

    delete_keys!(diagnostics_file, keys)
    save_diagnostics!(diagnostics_file, tracers, computed_output; eos)

    return nothing
end
function delete_keys!(diagnostics_file, keys)

    for k ∈ keys
        delete!(diagnostics_file, k)
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
function φ_interface_flux!(diagnostics_file::AbstractString, tracers::AbstractString, tracer::Symbol, eos)

    NCDataset(tracers) do ds

        φ = ds[tracer]
        Δφ₀ = φ[1, 1, 1, 1] - 0.5 * (φ[1, 1, 1, 1] - φ[1, 1, end, 1])
        timestamps = ds[:time][:]
        Δt = diff(timestamps)
        ΔV = diff(ds[:xC][1:2])[1] * diff(ds[:yC][1:2])[1] * diff(ds[:zC][1:2])[1]
        V = (1:length(reshape(φ[:, :, :, 1], :))) * ΔV
        SA = 0.07^2
        z✶ = V / SA
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

        save_fluxes!(diagnostics_file, φ_interface_flux, interface_idx, tracer, eos)

    end

    return nothing
end
"""
    function save_fluxes!(diagnostics_file, φ_interface_flux, interface_idx, tracer)
Save the flux through, and index of, the diffusive interface. **Note** the index is for the reshaped
and resorted vector.
"""
function save_fluxes!(diagnostics_file, φ_interface_flux, interface_idx, tracer, eos)

    if isfile(diagnostics_file)
        jldopen(diagnostics_file, "a+") do file
            file[eos*string(tracer)*"_flux"] = φ_interface_flux
            file[eos*string(tracer)*"_interface_idx"] = interface_idx
        end
    else
        jldopen(diagnostics_file, "w") do file
            file[eos*string(tracer)*"_flux"] = φ_interface_flux
            file[eos*string(tracer)*"_interface_idx"] = interface_idx
        end
    end

    return nothing
end
"""
    function interface_thickness!(diagnostics_file::AbstractString, tracers::AbstractString)
Calculate the interface thickness using the method in appendix A of
[Sommer et al. (2013)](https://journals.ametsoc.org/view/journals/atot/30/8/jtech-d-12-00272_1.xml).
"""
function interface_thickness!(diagnostics_file::AbstractString, tracers::AbstractString, eos)

    ds = NCDataset(tracers) do ds

        timestamps = ds[:time][:]
        ΔT = abs.(ds[:Tₗ_mean][:] .- ds[:Tᵤ_mean][:])
        ΔS = abs.(ds[:Sₗ_mean][:] .- ds[:Sᵤ_mean][:])
        zC = abs.(ds[:zC][:])

        Tmidpoint = 0.5 * (ds[:Tₗ_mean][1] .+ ds[:Tᵤ_mean][1])
        Smidpoint = 0.5 * (ds[:Sₗ_mean][1] .+ ds[:Sᵤ_mean][1])

        hₜ = similar(timestamps[1:end-1])
        hₛ = similar(timestamps[1:end-1])

        for t ∈ eachindex(timestamps)
            T = reshape(mean(ds[:T][:, :, :, t+1], dims = (1, 2)), :)
            find = findall(Tmidpoint - (ΔT[t+1]/4) .<  T .< Tmidpoint + (ΔT[t+1]/4))
            intercept, slope = [ones(length(find)) zC[find]] \ T[find]
            hₜ[t] = ΔT[t] / slope
        end

        for t ∈ eachindex(timestamps)
            S = reshape(mean(ds[:S][:, :, :, t+1], dims = (1, 2)), :)
            find = findall(Smidpoint - (ΔS[t+1]/4) .<  S .< Smidpoint + (ΔS[t+1]/4))
            intercept, slope = [ones(length(find)) zC[find]] \ S[find]
            hₛ[t] = ΔS[t] / slope
        end

        r = hₜ ./ hₛ
        save_interface_thickness!(diagnostics_file, hₜ, hₛ, r, ΔS, ΔT, eos)
    end

    return nothing
end
function save_interface_thickness!(diagnostics_file, hₜ, hₛ, r, ΔS, ΔT, eos)

    if isfile(diagnostics_file)
        jldopen(diagnostics_file, "a+") do file
            file[eos*"hₜ"] = hₜ
            file[eos*"hₛ"] = hₛ
            file[eos*"r"] = r
            file[eos*"ΔS"] = ΔS
            file[eos*"ΔT"] = ΔT
        end
    else
        jldopen(diagnostics_file, "w") do file
            file[eos*"hₜ"] = hₜ
            file[eos*"hₛ"] = hₛ
            file[eos*"r"] = r
            file[eos*"ΔS"] = ΔS
            file[eos*"ΔT"] = ΔT
        end
    end

    return nothing
end
"""
    function Ẽ(diagnostics_file::AbstractString, co::AbstractString, tracers::AbstractString, eos;
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
The kwarg `interface_offset` is how far away to move (vertically) from the interface to define
the upper and lower layers.
"""
function Ẽ(diagnostics_file::AbstractString, co::AbstractString, tracers::AbstractString, eos;
            interface_offset = 1)

    ds_co = NCDataset(co)
    ds_tracers = NCDataset(tracers)

    timestamps = ds_co[:time][:]
    zC = ds_co[:zC][:]
    Δx, Δy, Δz = diff(ds_co[:xC][:])[1], diff(ds_co[:yC][:])[1], diff(ds_co[:zC][:])[1]

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
        lower = 1:interface_idx-interface_offset
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

    save_Ẽ!(diagnostics_file, Ẽ, T̄_ts, S̄_ts, ρ̄_ts, eos)

    return nothing
end
function save_Ẽ!(diagnostics_file, Ẽ, T̄_ts, S̄_ts, ρ̄_ts, eos)

    if isfile(diagnostics_file)
        jldopen(diagnostics_file, "a+") do file
            file[eos*"Ẽ"] = Ẽ
            file[eos*"Tₗ_Tᵤ_ts"] = T̄_ts
            file[eos*"Sₗ_Sᵤ_ts"] = S̄_ts
            file[eos*"ρₗ_ρᵤ_ts"] = ρ̄_ts
        end
    else
        jldopen(diagnostics_file, "w") do file
            file[eos*"Ẽ"] = Ẽ
            file[eos*"Tₗ_Tᵤ_ts"] = T̄_ts
            file[eos*"Sₗ_Sᵤ_ts"] = S̄_ts
            file[eos*"ρₗ_ρᵤ_ts"] = ρ̄_ts
        end
    end
    return nothing
end
