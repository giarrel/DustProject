using CairoMakie

const GM = 1.327e20::Float64                            # G * M_sol [m^3 / s^2]
const AU = 1.496e11::Float64                            # 1AU [m]
const yr = 3.154e7::Float64                             # 1yr [s]
const day = 24.0 * 60.0 * 60.0::Float64                 # 1d [s]

# Parameter definieren
q = 1.0
B = 1e-9
m = 1.0
w = q*B/m

# Anfangsbedingungen
v0 = [ 30e3, 0.0,  30e3] # [vx(0), vy(0), vz(0)]
x0 = [0.0, 0.0, 0.0] # [x(0), y(0), z(0)]

# Zeitbereich
t = 0:day:600yr

#= Berechnung der Positionen
x = x0[1] .+ v0[1]/w*sin.(w*t) .- v0[2]/w*(1 .- cos.(w*t))
y = x0[2] .+ v0[2]/w*sin.(w*t) .+ v0[1]/w*(1 .- cos.(w*t))
z = x0[3] .+ v0[3]*t
=#
x = v0[1]/w*sin.(w*t)
y = v0[1]/w*(1 .- cos.(w*t))
z = v0[3]*t

# Positionen plotten
fig = Figure()
ax = Axis3(fig[1, 1]; xlabel = "AU", ylabel = "AU", zlabel = "AU", title = "Magnetic Field Trajectory Analytical Solution")
lines!(ax, Float32.(x)/AU, Float32.(y)/AU, Float32.(z)/AU, markersize = 0.1)
fig
save("DustProject/comet_plots/magnetic_trajectory.pdf", fig)