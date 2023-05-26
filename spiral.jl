using Plots

const Tsolarrotate = 25.05 * 24 * 3600  # [s]
const ωsun = 2π / Tsolarrotate  # [rad/s]
const Vr = 300e3  # [m/s]

function parker_spiral_2d(rmin, rmax, num_points=300)
    φ₀ = 0
    r = range(rmin, rmax, length=num_points)
    φ = @. 2π * (φ₀ + rmin * ωsun / Vr * (r / rmin - log(r) - 1 + log(rmin)))

    return r, φ
end
function plot_parker_spiral_2d()
    rmin = 1.0 * Rsun  # inner radius, [m]
    rmax = 5.3 * AU    # outer radius, [m]

    r, φ = parker_spiral_2d(rmin, rmax)

    plot(φ, r./AU, proj=:polar, label="Parker spiral", title="Interplanetary Magnetic Field In the Ecliptic Plane", lw=2)
    plot!(0:0.1:2π, ones(63), linestyle=:dash, proj=:polar, label="Earth")


    plot!(xlims=(0, 2π), ylims=(0, rmax/AU), size=(800, 800))
    xlabel!("φ")
    ylabel!("Distance (AU)")

    return current()
end

spiral_plot = plot_parker_spiral_2d()
display(spiral_plot)