using DifferentialEquations, LinearAlgebra, CairoMakie

# Constants
const GM = 1.327e20
const AU = 1.496e11
const day = 24.0 * 60.0 * 60.0
const yr = 3.154e7
beta = 0
qm = 1
B = 0.000005
B_field = [0,0,B]
tspan = (0.0, 10e4*yr)
stepsize = 0.01 * day

# Function to compute the analytical solution
v_anal(r) = qm * r * B/2 + sqrt((qm*B*r)^2 / 4 + GM / r)

# Initial conditions
initial_pos = [AU,0,0]
initial_vel = [0, v_anal(norm(initial_pos)), 0]
u0 = vcat(initial_pos, initial_vel)

# Define the ODE system
function orbit_ode!(du, u, p, t)
    beta, GM = p
    r = norm(u[1:3])
    du[1:3] .= u[4:6]
    du[4:6] .= -GM * u[1:3] / r^3 * (1 - beta) - qm * cross(u[4:6], B_field) 
end

method_string = [
#    Euler(),
#    Midpoint(),
#    ImplicitEuler(),
#    ImplicitMidpoint(),
#    Trapezoid(),
    RK4(),
]

# Initialize a dictionary to store errors, times, and computation times for each method
results_dict = Dict{String, Tuple{Vector{Float64}, Vector{Float64}, Float64}}()
methods_used = ""

# Iterate over the list of integration methods
for method in method_string
    local comp_time = @elapsed begin
        # Solve the ODE system
        local prob = ODEProblem(orbit_ode!, u0, tspan, (beta, GM))
        local sol = solve(prob, method, adaptive=false, dt=stepsize)
    end

    local r_numeric = [norm(u[1:3]) for u in sol.u]
    local v_num = [norm(u[4:6]) for u in sol.u]

    local times = sol.t

    local rel_err = [abs(v_num[i]-v_anal(r_numeric[i]))/v_anal(r_numeric[i]) for i in 1:length(times)]

    # Method name for labels
    local method_label = split(split(string(method), '{')[1], '(')[1]

    # Store the errors, times, and computation time for the current method
    global results_dict[method_label] = (rel_err, times, comp_time)

    # Add the name of the current method to the string
    global methods_used *= method_label * "_"
end

# Remove the trailing underscore from the methods_used string
methods_used = chop(methods_used)

# Figure
f = Figure()
Axis(f[1, 1];yscale=log10,title = "Grav_mag_err_stepsize=$(stepsize/day) days", xlabel="Time [y]", ylabel="Velocity error relative [%]")

for (method_label, (error, times, comp_time)) in results_dict
    # Filter out the values that are too small for the log plot
    local valid_indices = filter(i -> error[i] > 0, 1:length(error))

    # Use only valid indices for plotting
    local filtered_error = error[valid_indices]
    local filtered_times = times[valid_indices] /day/365

    scatter!(filtered_times, filtered_error, label="$(method_label) ($(comp_time) s)")
end

axislegend()

# Save the plot
save("gravmagplots/Grav_mag_err_yr$(tspan[2]/yr)_stepsize$(stepsize/day)_days_methods_$(methods_used).png", f)
GC.gc()