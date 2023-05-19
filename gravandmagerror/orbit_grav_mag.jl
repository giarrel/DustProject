using DifferentialEquations, LinearAlgebra, CairoMakie

# Constants
const GM = 1.327e20
const AU = 1.496e11
const day = 24.0 * 60.0 * 60.0
const yr = 3.154e7
beta = 0
qm = 1
B = 0.005
B_field = [0,0,B]
tspan = (0.0, 10yr)
stepsize = 0.01 * day

# Function to compute the analytical solution
v_anal(r) = qm * r * B/2 + sqrt((qm*B*r)^2 / 4 + GM / r)

# Initial conditions
initial_pos = [10AU,0,0]
initial_vel = [0, v_anal(norm(initial_pos)), 0]
u0 = vcat(initial_pos, initial_vel)

# Define the ODE system
function orbit_ode!(du, u, p, t)
    beta, GM = p
    r = norm(u[1:3])
    du[1:3] .= u[4:6]
    du[4:6] .= -GM * u[1:3] / r^3 * (1 - beta) - qm * cross(u[4:6], B_field) 
end

# Solve the ODE system
prob = ODEProblem(orbit_ode!, u0, tspan, (beta, GM))
sol = solve(prob, Trapezoid(), adaptive=false, dt=stepsize)

# Extract the position components
x_pos_numeric = [u[1] for u in sol.u]
y_pos_numeric = [u[2] for u in sol.u]
z_pos_numeric = [u[3] for u in sol.u]

x_vel_numeric = [point[4] for point in sol.u]
y_vel_numeric = [point[5] for point in sol.u]
z_vel_numeric = [point[6] for point in sol.u]

# Plot the results
figure = lines(x_pos_numeric/AU,y_pos_numeric/AU)
save("comet_plots/circular_orbit_magn_grav_B$(B)_R$(norm(initial_pos)/AU)_t$(tspan[2])).png", figure)
