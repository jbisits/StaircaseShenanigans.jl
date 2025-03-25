"""
    function StaircaseShenanigans.animate_tracers(tracers::AbstractString)
Animate the salinity and temperature `tracers` from saved `.nc` output.
"""
function StaircaseShenanigans.animate_tracers(tracers::AbstractString; xslice = 52, yslice = 52)

    NCDataset(tracers) do ds

        x = ds["x_caa"][:]
        z = ds["z_aac"][:]
        t = ds["time"][:]

        n = Observable(1)
        S = @lift ds["S"][:, yslice, :, $n]
        S_profile = @lift ds["S"][xslice, yslice, :, $n]
        Θ = @lift ds["T"][:, yslice, :, $n]
        Θ_profile = @lift ds["T"][xslice, yslice, :, $n]
        time_title = @lift @sprintf("t=%1.2f minutes", t[$n] / 60)

        fig = Figure(size = (1000, 1000))
        ax = [Axis(fig[j, i], title = (i == 1 && j == 1) ? time_title : "") for i ∈ 1:2, j ∈ 1:2]

        # Salinity
        lines!(ax[1], S_profile, z)
        ax[1].xlabel = "S gkg⁻¹"
        ax[1].ylabel = "z (m)"
        ax[1].xaxisposition = :top

        Scmap = cgrad(:haline)[2:end-1]
        Srange = extrema(ds[:S][:, :, :, 1])
        Slow = cgrad(:haline)[1]
        Shigh = cgrad(:haline)[end]
        hmS = heatmap!(ax[2], x, z, S, colorrange = Srange, colormap = Scmap,
                        lowclip = Slow, highclip = Shigh)

        ax[2].xlabel = "x (m)"
        ax[2].ylabel = "z (m)"
        Colorbar(fig[1, 3], hmS, label = "S gkg⁻¹")

        linkyaxes!(ax[1], ax[2])
        hideydecorations!(ax[2], ticks = false)

        # Temperature
        lines!(ax[3], Θ_profile, z)
        ax[3].xlabel = "Θ°C"
        ax[3].ylabel = "z (m)"
        ax[3].xaxisposition = :top

        Θcmap = cgrad(:thermal)[2:end-1]
        Θrange = extrema(ds[:T][:, :, :, 1])
        Θlow = cgrad(:thermal)[1]
        Θhigh = cgrad(:thermal)[end]
        hmΘ = heatmap!(ax[4], x, z, Θ, colorrange = Θrange, colormap = Θcmap,
                        lowclip = Θlow, highclip = Θhigh)

        ax[4].xlabel = "x (m)"
        ax[4].ylabel = "z (m)"
        Colorbar(fig[2, 3], hmΘ, label = "Θ°C")

        linkyaxes!(ax[3], ax[4])
        hideydecorations!(ax[4], ticks = false)

        frames = eachindex(t)
        record(fig, joinpath(pwd(), "tracers.mp4"),
            frames, framerate=8) do i
            msg = string("Plotting frame ", i, " of ", frames[end])
            print(msg * " \r")
            n[] = i
        end

    end

    return nothing
end
"""
    function StaircaseShenanigans.animate_tracers_anomaly(tracers::AbstractString)
Animate the salinity and temperature `tracers` anomalies from saved `.nc` output.
"""
function StaircaseShenanigans.animate_tracers_anomaly(tracers::AbstractString; xslice = 52, yslice = 52)

    NCDataset(tracers) do ds

        x = ds["x_caa"][:]
        z = ds["z_aac"][:]
        t = ds["time"][:]

        n = Observable(1)
        S = @lift ds[:S′][:, yslice, :, $n]
        S_profile = @lift ds[:S′][xslice, yslice, :, $n]
        Θ = @lift ds[:T′][:, yslice, :, $n]
        Θ_profile = @lift ds[:T′][xslice, yslice, :, $n]
        time_title = @lift @sprintf("t=%1.2f minutes", t[$n] / 60)

        fig = Figure(size = (1000, 1000))
        ax = [Axis(fig[j, i], title = (i == 1 && j == 1) ? time_title : "") for i ∈ 1:2, j ∈ 1:2]

        # Salinity
        Srange = extrema(ds[:S′][:, :, :, end])
        lines!(ax[1], S_profile, z)
        ax[1].xlabel = "S′ gkg⁻¹"
        ax[1].ylabel = "z (m)"
        ax[1].xaxisposition = :top
        xlims!(ax[1], Srange)

        Scmap = cgrad(:haline)[2:end-1]
        Slow = cgrad(:haline)[1]
        Shigh = cgrad(:haline)[end]
        hmS = heatmap!(ax[2], x, z, S, colorrange = Srange, colormap = Scmap,
                        lowclip = Slow, highclip = Shigh)

        ax[2].xlabel = "x (m)"
        ax[2].ylabel = "z (m)"
        Colorbar(fig[1, 3], hmS, label = "S′ gkg⁻¹")

        linkyaxes!(ax[1], ax[2])
        hideydecorations!(ax[2], ticks = false)

        # Temperature
        Θrange = extrema(ds[:T′][:, :, :, end])
        lines!(ax[3], Θ_profile, z)
        ax[3].xlabel = "Θ°C"
        ax[3].ylabel = "z (m)"
        ax[3].xaxisposition = :top
        xlims!(ax[3], Θrange)

        Θcmap = cgrad(:thermal)[2:end-1]
        Θlow = cgrad(:thermal)[1]
        Θhigh = cgrad(:thermal)[end]
        hmΘ = heatmap!(ax[4], x, z, Θ, colorrange = Θrange, colormap = Θcmap,
                        lowclip = Θlow, highclip = Θhigh)

        ax[4].xlabel = "x (m)"
        ax[4].ylabel = "z (m)"
        Colorbar(fig[2, 3], hmΘ, label = "Θ°C")

        linkyaxes!(ax[3], ax[4])
        hideydecorations!(ax[4], ticks = false)

        frames = eachindex(t)
        record(fig, joinpath(pwd(), "tracers_anomaly.mp4"),
            frames, framerate=8) do i
            msg = string("Plotting frame ", i, " of ", frames[end])
            print(msg * " \r")
            n[] = i
        end

    end

    return nothing
end
"""
    function animate_density(computed_output::AbstractString, variable::AbstractString=σ;
                                     xslice = 52, yslice = 52)
Animate the density variable in `computed_output`.
"""
function StaircaseShenanigans.animate_density(computed_output::AbstractString, variable::AbstractString;
                               xslice = 52, yslice = 52)

    NCDataset(computed_output) do ds

        x = ds["x_caa"][:]
        z = ds["z_aac"][:]
        t = ds["time"][:]

        n = Observable(1)
        σ = @lift ds[variable][:, yslice, :, $n]
        σ_profile = @lift ds[variable][xslice, yslice, :, $n]
        time_title = @lift @sprintf("t=%1.2f minutes", t[$n] / 60)

        fig = Figure(size = (1000, 500))
        ax = [Axis(fig[1, i], title = i == 1 ? time_title : "") for i ∈ 1:2]

        lines!(ax[1], σ_profile, z)
        ax[1].xlabel = "σ₀ kgm⁻³"
        ax[1].ylabel = "z"
        ax[1].xaxisposition = :top
        ax[1].xticklabelrotation = π / 4
        xlims!(ax[1], extrema(ds[variable][:, :, :, end]))

        colormap = cgrad(:dense)[2:end-1]
        colorrange = extrema(ds[variable][:, :, :, end])
        lowclip = cgrad(:dense)[1]
        highclip = cgrad(:dense)[end]
        hm = heatmap!(ax[2], x, z, σ; colorrange, colormap, lowclip, highclip)

        ax[2].xlabel = "x (m)"
        ax[2].ylabel = "z (m)"
        Colorbar(fig[1, 3], hm, label = "σ₀ kgm⁻³")

        linkyaxes!(ax[1], ax[2])
        hideydecorations!(ax[2], ticks = false)

        frames = eachindex(t)
        record(fig, joinpath(pwd(), "density.mp4"),
            frames, framerate=8) do i
            msg = string("Plotting frame ", i, " of ", frames[end])
            print(msg * " \r")
            n[] = i
        end

    end

    return nothing
end
"""
    function animate_density_anomaly(computed_output::AbstractString, variable::AbstractString=σ;
                                     xslice = 52, yslice = 52)
Animate the density anomaly from `variable` and the background `variable` in `computed_output`.
"""
function StaircaseShenanigans.animate_density_anomaly(computed_output::AbstractString, variable::AbstractString;
                               xslice = 52, yslice = 52)

    NCDataset(computed_output) do ds

        x = ds["x_caa"][:]
        z = ds["z_aac"][:]
        t = ds["time"][:]

        n = Observable(1)
        σ_backgroud = ds[variable*"_background"][:, yslice, :]
        σ = @lift ds[variable][:, yslice, :, $n] .- σ_backgroud
        σ_backgroud_profile =  ds[variable*"_background"][xslice, yslice, :]
        σ_profile = @lift ds[variable][xslice, yslice, :, $n] .- σ_backgroud_profile
        time_title = @lift @sprintf("t=%1.2f minutes", t[$n] / 60)

        fig = Figure(size = (1000, 500))
        ax = [Axis(fig[1, i], title = i == 1 ? time_title : "") for i ∈ 1:2]

        colorrange = extrema(ds[variable][:, :, :, end] .- ds[variable*"_background"][:, :, :])

        lines!(ax[1], σ_profile, z)
        ax[1].xlabel = "σ₀′ kgm⁻³"
        ax[1].ylabel = "z"
        ax[1].xaxisposition = :top
        ax[1].xticklabelrotation = π / 4
        xlims!(ax[1], colorrange)

        colormap = cgrad(:dense)[2:end-1]
        lowclip = cgrad(:dense)[1]
        highclip = cgrad(:dense)[end]
        hm = heatmap!(ax[2], x, z, σ; colorrange, colormap, lowclip, highclip)

        ax[2].xlabel = "x (m)"
        ax[2].ylabel = "z (m)"
        Colorbar(fig[1, 3], hm, label = "σ₀′ kgm⁻³")

        linkyaxes!(ax[1], ax[2])
        hideydecorations!(ax[2], ticks = false)

        frames = eachindex(t)
        record(fig, joinpath(pwd(), "density_anomaly.mp4"),
            frames, framerate=8) do i
            msg = string("Plotting frame ", i, " of ", frames[end])
            print(msg * " \r")
            n[] = i
        end

    end

    return nothing
end
"""
    function StaircaseShenanigans.animate_profile_in_S_Θ_space(tracers::AbstractString;
                                                                anomaly = false,
                                                                xslice = 52, yslice = 52)
Animate the `S` and `T` profiles at `xslice`, `yslice` in salinity-temperature space.
Setting `anomaly = true` will use the saved temperature anomaly.
"""
function StaircaseShenanigans.animate_profile_in_S_Θ_space(tracers::AbstractString;
                                                            anomaly = false,
                                                            xslice = 52, yslice = 52)

    NCDataset(tracers) do ds

        t = ds["time"][:]

        tracer_names = anomaly == false ? (S = :S, T = :T) : (S = :S′, T = :T′)

        n = Observable(1)
        S_profile = @lift ds[tracer_names.S][xslice, yslice, :, $n]
        T_profile = @lift ds[tracer_names.T][xslice, yslice, :, $n]
        time_title = @lift @sprintf("t=%1.2f minutes", t[$n] / 60)

        fig = Figure(size = (500, 500))
        ax = Axis(fig[1, 1], title = time_title, xlabel = "S (gkg⁻ꜝ)", ylabel = "Θ (°C)")

        scatter!(ax, S_profile, T_profile)

        frames = eachindex(t)
        record(fig, joinpath(pwd(), "tracerprofiles_in_S_T_space.mp4"),
            frames, framerate=8) do i
            msg = string("Plotting frame ", i, " of ", frames[end])
            print(msg * " \r")
            n[] = i
        end

    end

    return nothing
end
"""
    function visualise_initial_conditions(sdns::StaircaseDNS, xslice::Integer, yslice::Integer)
Plot the initial state of the `tracers` in a `model`. This function assumes there are two
tracers (salinity and temperature) and plots the x-z, y-z and field-z initial fields at
`xslice` and `yslice`.
"""
function StaircaseShenanigans.visualise_initial_conditions(sdns::StaircaseDNS, xslice::Integer, yslice::Integer)

    model = sdns.model

    x, y, z = nodes(model.grid, (Center(), Center(), Center()))
    S = model.tracers.S
    T = model.tracers.T
    fig = Figure(size = (600, 1500))
    ax = [Axis(fig[j, i]) for i ∈ 1:2, j ∈ 1:3]

    hm = heatmap!(ax[1], x, z, interior(S, :, yslice, :, 1); colormap = :haline)
    ax[1].title = "Initial salinity (x-z)"
    ax[1].xlabel = "x (m)"
    ax[1].ylabel = "z (m)"
    heatmap!(ax[2], x, z, interior(S, xslice, :, :, 1); colormap = :haline)
    ax[2].title = "Initial salinity (y-z)"
    ax[2].xlabel = "y (m)"
    ax[2].ylabel = "z (m)"
    Colorbar(fig[1, 3], hm, label = "S (gkg⁻¹)")
    hm = heatmap!(ax[3], x, z, interior(T, :, yslice, :, 1); colormap = :thermal)
    ax[3].title = "Initial temperature (x-z)"
    ax[3].xlabel = "x (m)"
    ax[3].ylabel = "z (m)"
    heatmap!(ax[4], y, z, interior(T, xslice, :, :, 1); colormap = :thermal)
    ax[4].title = "Initial temperature (y-z)"
    ax[4].xlabel = "y (m)"
    ax[4].ylabel = "z (m)"
    Colorbar(fig[2, 3], hm, label = "Θ (°C)")
    lines!(ax[5], interior(S, xslice, yslice, :, 1), z)
    ax[5].title = "Initial salinity profile"
    ax[5].xlabel = "S (gkg⁻¹)"
    ax[5].ylabel = "z (m)"
    lines!(ax[6], interior(T, xslice, yslice, :, 1), z)
    ax[6].title = "Initial temperature profile"
    ax[6].xlabel =  "Θ (°C)"
    ax[6].ylabel = "z (m)"

    return fig

end
"""
    function visualise_initial_density(sdns::StaircaseDNS, xslice::Integer,  yslice::Integer,
                                       pressure::Union{Number, Vector{Number}})
Compute and plot the initial density at `pressure` (either reference pressure or in-situ
pressure). The arguments `xslice` and `yslice` are used to choose where in the domain the
figures are from.
"""
function StaircaseShenanigans.visualise_initial_density(sdns::StaircaseDNS, xslice::Integer,  yslice::Integer,
                                                        geopotential_height::Integer)

    model = sdns.model

    x = xnodes(model.grid, Center(), Center(), Center())
    z = znodes(model.grid, Center(), Center(), Center())
    σ = Field(seawater_density(model; geopotential_height))
    compute!(σ)

    σ_hm = interior(σ, :, yslice, :)
    fig = Figure(size = (1000, 600))
    ax = [Axis(fig[1, i]) for i ∈ 1:2]

    hm = heatmap!(ax[1], x, z, σ_hm; colormap = :dense)
    ax[1].title = "Initial density (x-z)"
    ax[1].xlabel = "x (m)"
    ax[1].ylabel = "z (m)"
    Colorbar(fig[2, 1], hm, label = "ρ (kgm⁻³)", vertical = false, flipaxis = false)

    σ_profile = interior(σ, xslice, yslice, :)
    lines!(ax[2], σ_profile, z)
    ax[2].title = "Initial density profile"
    ax[2].xlabel = "ρ (kgm⁻³)"
    ax[2].ylabel = "z (m)"

    linkyaxes!(ax[1], ax[2])

    return fig

end
