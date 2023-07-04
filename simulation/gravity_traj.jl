using DifferentialEquations
using LinearAlgebra
using CairoMakie

include("constants.jl")
include("duststart.jl")

#paramerters
tspan=30yr
beta=0.3
no_of_traj = 5

function trajectory_ode(du, u, p, t)
        
    beta, GM = p
    r = norm(u[1:3])
    du[1:3] = u[4:6]
    du[4:6] = -GM * u[1:3] / r^3 * (1 - beta)

end

# Initial conditions and problem setup
p = (beta, GM)
u0 = startparams(no_of_traj)

# Definieren Sie eine Konditionsfunktion, die wahr wird, wenn die Trajektorie näher als 2.5 Sonnenradien an der Sonne ist.
function condition(dist,t,integrator) 
    s = integrator.u
    dist = sqrt(s[1]^2 + s[2]^2 + s[3]^2) # Die Entfernung vom Ursprung 
    dist < 3 *696340e3 # Die Bedingung für das Stoppen end
end

# Definieren Sie eine Funktionsfunktion, die das Abbruchereignis handhabt und den Integrator stoppt.
function affect!(integrator)
    terminate!(integrator)
end

# Erstellen Sie das Abbruchereignis mit den obigen Funktionen.
cb = DiscreteCallback(condition,affect!) # Erstellen Sie den Callback

fig = Figure()
ax = Axis3(fig[1, 1], perspectiveness = 0.4)

for i in axes(u0)[2]
    local prob = ODEProblem(trajectory_ode, u0[1:6,i], tspan, p, callback=cb)
    local sol = solve(prob, Vern9(), adaptive=false, dt=0.01day)

    # Extract the results
    local x_pos_numeric = [point[1] for point in sol.u]
    local y_pos_numeric = [point[2] for point in sol.u]
    local z_pos_numeric = [point[3] for point in sol.u]

    local x_vel_numeric = [point[4] for point in sol.u]
    local y_vel_numeric = [point[5] for point in sol.u]
    local z_vel_numeric = [point[6] for point in sol.u]

    local times = sol.t

    
    local energy_initial = 1/2*norm([x_vel_numeric[1],y_vel_numeric[1],z_vel_numeric[1]])^2-GM*(1 - beta)/norm([x_pos_numeric[1],y_pos_numeric[1],z_pos_numeric[1]])
    local energy_errors = [(1/2*norm([x_vel_numeric[j],y_vel_numeric[j],z_vel_numeric[j]])^2-GM*(1 - beta)/norm([x_pos_numeric[j],y_pos_numeric[j],z_pos_numeric[j]]) - energy_initial)/energy_initial for j in 1:length(times)]

    println("Maximum energy change ", maximum(abs.(energy_errors))*100, "%")
    

    lines!(ax, x_pos_numeric/AU, y_pos_numeric/AU, z_pos_numeric/AU, color = :blue, linewidth = 1)
end

scatter!(ax, [0], [0], [0], color = :orange, markersize = 10)
limits!(ax,-80,80,-80,80,-80,80)

ax.xlabel = "AU"
ax.ylabel = "AU"
ax.zlabel = "AU"


fig