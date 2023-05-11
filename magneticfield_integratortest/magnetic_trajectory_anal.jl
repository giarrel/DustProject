using GLMakie

# Parameter definieren
q = 1.0
B = 1.0
m = 1.0
w = q*B/m

# Anfangsbedingungen
v0 = [1.0, 0.0, 1.0] # [vx(0), vy(0), vz(0)]
x0 = [0.0, 0.0, 0.0] # [x(0), y(0), z(0)]

# Zeitbereich
t = 0:0.01:10

# Berechnung der Positionen
x = x0[1] .+ v0[1]/w*sin.(w*t) .- v0[2]/w*(1 .- cos.(w*t))
y = x0[2] .+ v0[2]/w*sin.(w*t) .+ v0[1]/w*(1 .- cos.(w*t))
z = x0[3] .+ v0[3]*t

# Positionen plotten
fig = Figure()
ax = Axis3(fig[1, 1]; title = "analy magn traj")
lines!(ax, Float32.(x), Float32.(y), Float32.(z), markersize = 0.1)
show(fig)