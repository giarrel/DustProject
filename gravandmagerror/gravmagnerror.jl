using DifferentialEquations, LinearAlgebra, CairoMakie

# Constants
const GM = 1.327e20
const AU = 1.496e11
const day = 24.0 * 60.0 * 60.0
const yr = 3.154e7
beta = 0
qm = 1
B = 1e-5
B_field = [0,0,B]
tspan = (0.0, 200yr)
stepsize = 0.1 * day
reltol=1e-10
adapt=true

# Function to compute the analytical solution
v_anal(r) = qm * r * B/2 + sqrt((qm*B*r)^2 / 4 + GM / r)
pos_anal(t,r) = r*[cos(v_anal(r)/r*t),sin(v_anal(r)/r*t),0]

# Initial conditions
initial_pos = [AU,0,0]
radius = norm(initial_pos)
initial_vel = [0, v_anal(radius), 0]
u0 = vcat(initial_pos, initial_vel)

# Define the ODE system
function orbit_ode!(du, u, p, t)
    beta, GM = p
    r = norm(u[1:3])
    du[1:3] .= u[4:6]
    du[4:6] .= -GM * u[1:3] / r^3 * (1 - beta) - qm * cross(u[4:6], B_field) 
end

method_string = [
    Euler(),

Midpoint(),
ImplicitEuler(),
ImplicitMidpoint(),
Trapezoid(),
RK4(),
Vern9(),
AutoVern9(Rodas5P())
    ]

# Initialize a dictionary to store errors, times, and computation times for each method
results_dict = Dict{String, Tuple{Vector{Float64}, Vector{Float64}, Float64}}()
methods_used = ""

# Iterate over the list of integration methods
for method in method_string
    local comp_time = @elapsed begin
        # Solve the ODE system
        local prob = ODEProblem(orbit_ode!, u0, tspan, (beta, GM))
        local sol = solve(prob, method, adaptive=adapt, reltol=reltol,dt=stepsize)
    end

    local pos_numeric = [u[1:3] for u in sol.u]

    local times = sol.t

    local rel_err_x = [abs((pos_anal(times[i],radius)[1] - pos_numeric[i][1])/pos_anal(times[i],radius)[1]) for i in 1:length(times)]
    local rel_err_y = [abs((pos_anal(times[i],radius)[2] - pos_numeric[i][2])/pos_anal(times[i],radius)[2]) for i in 1:length(times)]
    local rel_err_z = [abs(pos_anal(times[i],radius)[3] - pos_numeric[i][3]) / (abs(pos_anal(times[i],radius)[3]) + eps()) for i in 1:length(times)] #eps um zu vermeiden durch 0 teilen
    local rel_err = [norm([rel_err_x[i], rel_err_y[i], rel_err_z[i]]) for i in 1:length(times)]


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
Axis(f[1, 1];yscale=log10,title = "gravity and magnetic error, step size = $(stepsize/day) days", xlabel="Time [y]", ylabel="position error relative [%]")

for (method_label, (error, times, comp_time)) in results_dict
    # Filter out the values that are too small for the log plot
    local valid_indices = filter(i -> (Inf > error[i] > 0), 1:length(error))

    # Use only valid indices for plotting
    local filtered_error = error[valid_indices]
    local filtered_times = times[valid_indices] /day/365

    lines!(filtered_times, filtered_error, label="$(method_label) ($(comp_time) s)")
end

axislegend()

# Save the plot
#save("DustProject/error_plots_ada2/Grav_mag_(B=$(B)T)_(R=$(norm(initial_pos)/AU)AU)_integrationtime$(tspan[2]/yr)yr_reltol$(reltol)days_methods$(methods_used).pdf", f)
save("DustProject/error_plots_ada/Grav_mag_(B=$(B)T)_(R=$(norm(initial_pos)/AU)AU)_integrationtime$(tspan[2]/yr)yr_stepsize$(stepsize/day)days_methods$(methods_used).pdf", f)
GC.gc()