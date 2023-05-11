using DifferentialEquations, GLMakie

# Parameter definieren
q = 1.0
B = 1.0
m = 1.0
w = q*B/m

# Anfangsbedingungen
v0 = [1.0, 0.0, 1.0] # [vx(0), vy(0), vz(0)]
x0 = [0.0, 0.0, 0.0] # [x(0), y(0), z(0)]
u0 = vcat(x0, v0)

# Magnetfeld
B_field = [0, 0, B]

# Definition der Differentialgleichungen
function charged_particle!(du, u, p, t)
    du[1:3] = u[4:6] # dx/dt = v
    du[4] = q/m * (u[5]*B_field[3] - u[6]*B_field[2]) # dvx/dt = q/m * (vy*Bz - vz*By)
    du[5] = q/m * (u[6]*B_field[1] - u[4]*B_field[3]) # dvy/dt = q/m * (vz*Bx - vx*Bz)
    du[6] = q/m * (u[4]*B_field[2] - u[5]*B_field[1]) # dvz/dt = q/m * (vx*By - vy*Bx)
end

# Zeitbereich
tspan = (0.0, 10.0)

# Problem und Lösung
prob = ODEProblem(charged_particle!, u0, tspan)
sol = solve(prob)

# Positionen plotten
fig = Figure()
ax = Axis3(fig[1, 1], title = "Charged Particle Trajectory")
lines!(ax, Float32.(sol[1,:]), Float32.(sol[2,:]), Float32.(sol[3,:]), linewidth = 1, color = :blue)
fig
