using CairoMakie
using LinearAlgebra

include("KeptoCartfkt.jl")
include("cometlaunch.jl")
include("comnames.jl")
include("constants.jl")

#method
method_string = [
#KenCarp58(),
#RK4(),
Tsit5(),
Vern9(),
#VCABM5(),
#RadauIIA5(),
#Rodas5P()
#VCAB3(),
#VCAB4(),
#VCAB5(),
#VCABM3(),
#VCABM4(),
#VCABM5(),
#VCABM(),
#AN5(),
#JVODE_Adams()
]

beta=0
abstol=1e-9
reltol=1e-9

comet_name = "Phaethon" #Phaethon , Arend , Tuttle
# Initial conditions and problem setup
tspan = (0, 20000comets[comet_name].period)
initial_pos = keplerian_to_cartesian(comet_name, tspan[1], tspan[1])[1]
initial_vel = keplerian_to_cartesian(comet_name, tspan[1], tspan[1])[2]

# Initialize a dictionary to store errors, times, and computation times for each method
results_dict = Dict{String, Tuple{Vector{Float64}, Vector{Float64}, Float64}}()

# Iterate over the list of integration methods
for method in method_string
    local comp_time = @elapsed begin
        local vel_num, pos_num, times = launch_from_comet_1ord(initial_vel, initial_pos, tspan, beta, method,adaptive=true, abstol=abstol, reltol=reltol)
    end

    # Generate orbit points using the current times vector
    local positions = [keplerian_to_cartesian(comet_name, times[1], time)[1] for time in times]

    # Generate orbit points and compute errors
    local error_x = [abs((positions[i][1] - pos_num[1][i])/positions[i][1]) for i in 1:length(times)]
    local error_y = [abs((positions[i][2] - pos_num[2][i])/positions[i][2]) for i in 1:length(times)]
    local error_z = [abs((positions[i][3] - pos_num[3][i])/positions[i][3]) for i in 1:length(times)]

    local error = [norm([error_x[i], error_y[i], error_z[i]]) for i in 1:length(times)]
    local method_label = split(split(string(method), '{')[1], '(')[1]

    # Store the errors, times, and computation time for the current method
    results_dict[method_label] = (error, times, comp_time)
end

# Create the plot
f = Figure()
Axis(f[1, 1];yscale=log10,title = "1nd Order tol=$(reltol)", xlabel="Time [days]", ylabel="Error relative [%]")


for (method_label, (error, times, comp_time)) in results_dict
    # Filter out the values that are too small for the log plot
    valid_indices = filter(i -> error[i] > 0, 1:length(error))

    # Use only valid indices for plotting
    filtered_error = error[valid_indices]
    filtered_times = times[valid_indices] / day

    lines!(filtered_times, filtered_error, label="$(method_label) ($(comp_time) s)")
end

axislegend()

f