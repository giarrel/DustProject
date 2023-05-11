using CairoMakie
using LinearAlgebra

include("KeptoCartfkt.jl")
include("cometlaunch.jl")
include("comnames.jl")
include("constants.jl")

#method
method_string = [
#    SymplecticEuler(),
#    VelocityVerlet(),
#    VerletLeapfrog(),
    SofSpa10()
]

beta=0
stepsize=0.1day
periods=10000
last_N_timesteps_outside = 50000 #put 0 for all

comet_name = "Phaethon" #Phaethon , Arend , Tuttle

# Initial conditions and problem setup
tspan = (0, periods*comets[comet_name].period)
initial_pos = keplerian_to_cartesian(comet_name, tspan[1], tspan[1])[1]
initial_vel = keplerian_to_cartesian(comet_name, tspan[1], tspan[1])[2]

# Initialize a dictionary to store errors, times, and computation times for each method
results_dict = Dict{String, Tuple{Vector{Float64}, Vector{Float64}, Float64, Float64}}()

methods_used = ""

# Iterate over the list of integration methods
for method in method_string
    local comp_time = @elapsed begin
        local vel_num, pos_num, times = launch_from_comet(initial_vel, initial_pos, tspan, beta, method, stepsize=stepsize)
    end

    local last_N_timesteps = last_N_timesteps_outside

    #methoden die abbrechen funktionieren nicht mit nlasttimesteps logik
    if last_N_timesteps != 0
        if length(times) < (tspan[2]/stepsize)
            println("skipped")
            continue
        end
    end

    println("not skipped")

    # Check if last_N_timesteps is 0, if so, use the full length of times
    last_N_timesteps = last_N_timesteps == 0 ? length(times) : last_N_timesteps

    # Generate orbit points using the current times vector
    local comp_time2 = @elapsed begin
        local positions = [keplerian_to_cartesian(comet_name, times[1], time)[1] for time in times[end-last_N_timesteps+1:end]]
    end
    # Generate orbit points and compute errors
    local error_x = [abs((positions[i][1] - pos_num[1][i+length(times)-last_N_timesteps])/positions[i][1]) for i in 1:last_N_timesteps]
    local error_y = [abs((positions[i][2] - pos_num[2][i+length(times)-last_N_timesteps])/positions[i][2]) for i in 1:last_N_timesteps]
    local error_z = [abs((positions[i][3] - pos_num[3][i+length(times)-last_N_timesteps])/positions[i][3]) for i in 1:last_N_timesteps]

    local error = [norm([error_x[i], error_y[i], error_z[i]]) for i in 1:last_N_timesteps]
    local method_label = split(split(string(method), '{')[1], '(')[1]

    # Store the errors, times, and computation time for the current method
    results_dict[method_label] = (error, times[end-last_N_timesteps+1:end], comp_time, comp_time2)

    # Add the name of the current method to the string
    global methods_used *= method_label * "_"
end

methods_used = chop(methods_used)

# Create the plot
f = Figure()
Axis(f[1, 1];yscale=log10,title = "symplekt stepsize=$(stepsize/day) days", xlabel="Time [y]", ylabel="Error relative [%]")


for (method_label, (error, times, comp_time, comp_time2)) in results_dict
    # Filter out the values that are too small for the log plot
    valid_indices = filter(i -> error[i] > 0, 1:length(error))

    # Use only valid indices for plotting
    filtered_error = error[valid_indices]
    filtered_times = times[valid_indices] / day/365

    lines!(filtered_times, filtered_error, label="$(method_label) ($(comp_time) s)")
end

axislegend()

# Save the plot
save("error_plots/sympl_$(comet_name)_stepsize$(stepsize/day)_days_periods$(periods)_methods_$(methods_used).png", f)
GC.gc()