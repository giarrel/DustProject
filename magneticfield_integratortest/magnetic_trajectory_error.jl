using DifferentialEquations, CairoMakie, LinearAlgebra

#constants
const GM = 1.327e20::Float64                            # G * M_sol [m^3 / s^2]
const AU = 1.496e11::Float64                            # 1AU [m]
const yr = 3.154e7::Float64                             # 1yr [s]
const day = 24.0 * 60.0 * 60.0::Float64                 # 1d [s]

# Parameter definieren
q = 1.0
B = 5e-9
m = 1.0
w = q*B/m

# Anfangsbedingungen
v0 = [30e3, 0.0, 30e3] # [vx(0), vy(0), vz(0)]
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
tspan = (0.0, 1000yr)

# Problem und Lösung
prob = ODEProblem(charged_particle!, u0, tspan)
sol = solve(prob,RK4(),dt=day,adaptive=false)

# Analytische Lösung
x_analytic(t) = x0[1] .+ v0[1]/w*sin.(w*t) .- v0[2]/w*(1 .- cos.(w*t))
y_analytic(t) = x0[2] .+ v0[2]/w*sin.(w*t) .+ v0[1]/w*(1 .- cos.(w*t))
z_analytic(t) = x0[3] .+ v0[3]*t

# Fehler berechnen
error_x = abs.((x_analytic(sol.t) .- sol[1,:]) ./ x_analytic(sol.t))
error_y = abs.((y_analytic(sol.t) .- sol[2,:]) ./ y_analytic(sol.t))
error_z = abs.((z_analytic(sol.t) .- sol[3,:]) ./ z_analytic(sol.t))

# Betrag des relativen Fehlervektors
error_magnitude = sqrt.(error_x.^2 .+ error_y.^2 .+ error_z.^2)

# Fehler für jeden Zeitschritt plotten
fig = Figure()
ax = Axis(fig[1, 1],title = "magnetic error", xlabel="Time [y]", ylabel="Error relative [%]")
lines!(ax, sol.t/yr, error_magnitude)
fig