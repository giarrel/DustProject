using CairoMakie
using LinearAlgebra

include("KeptoCartfkt.jl")
include("cometlaunch.jl")
include("comnames.jl")
include("constants.jl")

#method
method_string = [
    

Euler(),
Midpoint(),
ImplicitEuler(),
ImplicitMidpoint(),
Trapezoid(),
RK4(),
Vern9(),
AutoVern9(Rodas5P())
#=
    AB3(),
    AB4(),
    AB5(),
    ABM32(),
    ABM43(),
    ABM54()
    =#
]

beta=0

reltol=1e-10
periods=150
last_N_timesteps_outside = 000 #0 for all

comet_name = "Phaethon" #Phaethon , Arend , Tuttle
# Initial conditions and problem setup
tspan = (0, periods*comets[comet_name].period)
initial_pos = keplerian_to_cartesian(comet_name, tspan[1], tspan[1])[1]
initial_vel = keplerian_to_cartesian(comet_name, tspan[1], tspan[1])[2]

# Initialize a dictionary to store errors, times, and computation times for each method
results_dict = Dict{String, Tuple{Vector{Float64}, Vector{Float64}, Float64, Float64}}()

# Create a string to store the names of the used methods
methods_used = ""

# Iterate over the list of integration methods
for method in method_string
    local comp_time = @elapsed begin
        local vel_num, pos_num, times = launch_from_comet_1ord(initial_vel, initial_pos, tspan, beta, method, adaptive=true, reltol=reltol)#abstol=abstol,
    end

    local last_N_timesteps = last_N_timesteps_outside #damit nicht vom ersten

    #methoden die abbrechen funktionieren nicht mit nlasttimesteps logik
    if last_N_timesteps != 0
        if times[end] < tspan[end]
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

    #= Generate orbit points and compute errors
    local error_x = [abs((positions[i][1] - pos_num[1][i])/positions[i][1]) for i in 1:last_N_timesteps]
    local error_y = [abs((positions[i][2] - pos_num[2][i])/positions[i][2]) for i in 1:last_N_timesteps]
    local error_z = [abs((positions[i][3] - pos_num[3][i])/positions[i][3]) for i in 1:last_N_timesteps]
=#
    local error_x = [abs((positions[i][1] - pos_num[1][i])) for i in 1:last_N_timesteps]
    local error_y = [abs((positions[i][2] - pos_num[2][i])) for i in 1:last_N_timesteps]
    local error_z = [abs((positions[i][3] - pos_num[3][i])) for i in 1:last_N_timesteps]

    local error = [norm([error_x[i], error_y[i], error_z[i]]) for i in 1:last_N_timesteps]
    local method_label = split(split(string(method), '{')[1], '(')[1]

    # Store the errors, times, and computation time for the current method
    results_dict[method_label] = (error, times[end-last_N_timesteps+1:end], comp_time, comp_time2)

    # Add the name of the current method to the string
    global methods_used *= method_label * "_"
end

# Remove the trailing underscore from the methods_used string
methods_used = chop(methods_used)

# Create the plot
f = Figure()
Axis(f[1, 1];yscale=log10,title = "Phaethon's orbit error, adaptive steps with tolerance = $(reltol)", xlabel="Time [y]", ylabel="Error relative [%]")


for (method_label, (error, times, comp_time, comp_time2)) in results_dict
    # Filter out the values that are too small for the log plot
    valid_indices = filter(i -> error[i] > 0, 1:length(error))

    # Use only valid indices for plotting
    filtered_error = error[valid_indices]
    filtered_times = times[valid_indices] /day/365

    lines!(filtered_times, filtered_error, label="$(method_label) ($(comp_time) s)")
end

axislegend()

# Save the plot
save("DustProject/error_plots_ada_abserror/1ode_$(comet_name)_tol$(reltol)_days_periods$(periods)_methods_$(methods_used).pdf", f)
GC.gc()