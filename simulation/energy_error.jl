using DifferentialEquations
using LinearAlgebra
using Plots

include("constants.jl")
include("duststart.jl")

#parameters
tspan = 30yr
beta = 0
no_of_traj = 10

function trajectory_ode(du, u, para, t)
    beta, GM = para
    r = norm(u[1:3])
    du[1:3] = u[4:6]
    du[4:6] = -GM * u[1:3] / r^3 * (1 - beta)
end

# Initial conditions and problem setup
para = (beta, GM)
u0 = vcat(startparams(no_of_traj))

# Define a condition function that becomes true when the trajectory is closer than 2.5 solar radii to the sun.
function condition(dist,t,integrator) 
    s = integrator.u
    dist = sqrt(s[1]^2 + s[2]^2 + s[3]^2) # The distance from the origin 
    dist < 3*696340e3 # The condition for stopping
end

# Define a function function that handles the termination event and stops the integrator.
function affect!(integrator)
    terminate!(integrator)
end

# Create the termination event with the above functions.
cb = DiscreteCallback(condition,affect!)

# Create an empty plot
p = plot()

for i in axes(u0)[2]
    local prob = ODEProblem(trajectory_ode, u0[1:6,i], tspan, para, callback=cb)
    local sol = solve(prob, RK4(), adaptive=false, dt=0.01day)

    # Extract the results
    local x_pos_numeric = [point[1] for point in sol.u]
    local y_pos_numeric = [point[2] for point in sol.u]
    local z_pos_numeric = [point[3] for point in sol.u]

    local x_vel_numeric = [point[4] for point in sol.u]
    local y_vel_numeric = [point[5] for point in sol.u]
    local z_vel_numeric = [point[6] for point in sol.u]

    local times = sol.t

    # Define energy errors and distances
    local energy_initial = 1/2*norm([x_vel_numeric[1],y_vel_numeric[1],z_vel_numeric[1]])^2-GM*(1 - beta)/norm([x_pos_numeric[1],y_pos_numeric[1],z_pos_numeric[1]])
    local energy_errors = [(1/2*norm([x_vel_numeric[j],y_vel_numeric[j],z_vel_numeric[j]])^2-GM*(1 - beta)/norm([x_pos_numeric[j],y_pos_numeric[j],z_pos_numeric[j]]) - energy_initial)/energy_initial for j in 1:length(times)]

    # Add a small value to avoid taking log of 0
    energy_errors .+= 1e-20

    println("Maximum energy change ", maximum(abs.(energy_errors))*100, "%")

    # Plot the energy change error against distance to the sun for each trajectory
    local distances = [norm([x_pos_numeric[j],y_pos_numeric[j],z_pos_numeric[j]]) for j in 1:length(times)]

    # Take the log of absolute error and distance
    energy_errors_log = log10.(abs.(energy_errors))
    distances_log = log10.(distances)

    # Plot
    plot!(p, distances_log, energy_errors_log, legend = false)
end

title!(p, "Log-Log plot of relative energy change error")
xlabel!(p, "Log of distance to sun (m)")
ylabel!(p, "Log of absolute energy change error")

# Show the plot
p
