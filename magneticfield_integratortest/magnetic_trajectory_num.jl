using DifferentialEquations, GLMakie, LinearAlgebra

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

function lorenzf(q,m,v,B_field)
    return -q/m*cross(v,B_field)
end

# Definition der Differentialgleichungen
function charged_particle!(du, u, p, t)
    du[1:3] = u[4:6]
    du[4:6] = lorenzf(q,m,u[4:6],B_field)
end

# Zeitbereich
tspan = (0.0, 10.0)

# Problem und Lösung
prob = ODEProblem(charged_particle!, u0, tspan)
sol = solve(prob,RK4(),dt=0.1,adaptive=false)

# Positionen plotten
fig = Figure()
ax = Axis3(fig[1, 1], title = "Charged Particle Trajectory")
lines!(ax, Float32.(sol[1,:]), Float32.(sol[2,:]), Float32.(sol[3,:]), linewidth = 1, color = :blue)
fig
