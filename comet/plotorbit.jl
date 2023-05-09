using CairoMakie

include("KeptoCartfkt.jl")
include("comnames.jl")
include("constants.jl")

function plot_comet_orbit(comet_name)
    local period = comets[comet_name].period
    local times_all = range(0, stop=1 * period, step=1 * day)

    # Generate orbit points for all times
    local positions_all = [keplerian_to_cartesian(comet_name, times_all[1], time)[1] for time in times_all]

    # Extract X, Y, and Z components of positions and convert to AU
    local x_positions = [pos[1] / AU for pos in positions_all]
    local y_positions = [pos[2] / AU for pos in positions_all]
    local z_positions = [pos[3] / AU for pos in positions_all]

    # Create the figure and axis
    fig = Figure(resolution = (800, 600))
    ax = Axis3(fig[1, 1], xlabel="X (AU)", ylabel="Y (AU)", zlabel="Z (AU)")

    # Plot the comet's orbit
    lines!(ax, x_positions, y_positions, z_positions, color = :blue, linewidth = 1.5, label = comet_name)

    # Plot the Sun at the origin as a red star
    scatter!(ax, [0], [0], [0], color = :red, markersize = 10, markershape = :star5, label = "Sun")

    # Add the legend
    axislegend(ax)

    # Save the plot
    save("comet_plots/$(comet_name)_orbit.png", fig)

    println("Plot saved as 'comet_plots/$(comet_name)_orbit.png'")
end

# Iterate over all comets and plot their orbits
for comet_name in keys(comets)
    plot_comet_orbit(comet_name)
end

