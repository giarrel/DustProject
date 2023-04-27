using Plots
using LinearAlgebra

include("KeptoCartfkt.jl")
include("cometlaunch.jl")
include("comnames.jl")

#constants
const GM = 1.327e20::Float64                            # G * M_sol [m^3 / s^2]
const AU = 1.496e11::Float64                            # 1AU [m]
const yr = 3.154e7::Float64                             # 1yr [s]
const day = 24.0 * 60.0 * 60.0::Float64                 # 1d [s]

#method
method_string = [
#    SymplecticEuler(),
#    VelocityVerlet(),
#    VerletLeapfrog(),
#    PseudoVerletLeapfrog(),
#    McAte2(),
#    Ruth3(),
#    McAte3(),
#    CandyRoz4(),
#    McAte4(),
    CalvoSanz4(),
#    McAte42(),
#    McAte5(),
    Yoshida6(),
    KahanLi6(),
    McAte8(),
    KahanLi8(),
    SofSpa10()
]

beta=0

comet_name = "Phaethon"
# Initial conditions and problem setup
tspan = (0, 2comets[comet_name].period)
initial_pos = keplerian_to_cartesian(comet_name, tspan[1], tspan[1])[1]
initial_vel = keplerian_to_cartesian(comet_name, tspan[1], tspan[1])[2]

# Initialize a dictionary to store errors, times, and computation times for each method
results_dict = Dict{String, Tuple{Vector{Float64}, Vector{Float64}, Float64}}()

# Iterate over the list of integration methods
for method in method_string
    local comp_time = @elapsed begin
        local vel_num, pos_num, times = launch_from_comet(initial_vel, initial_pos, tspan, beta, method)
    end

    # Generate orbit points using the current times vector
    local positions = [keplerian_to_cartesian(comet_name, times[1], time)[1] for time in times]

    # Generate orbit points and compute errors
    local error_x = [abs(positions[i][1] - pos_num[1][i]) for i in 1:length(times)]
    local error_y = [abs(positions[i][2] - pos_num[2][i]) for i in 1:length(times)]
    local error_z = [abs(positions[i][3] - pos_num[3][i]) for i in 1:length(times)]

    local error = [norm([error_x[i], error_y[i], error_z[i]]) for i in 1:length(times)]
    local method_label = split(split(string(method), '{')[1], '(')[1]

    # Store the errors, times, and computation time for the current method
    results_dict[method_label] = (error, times, comp_time)
end

# Create the plot
p = plot(title="symplectic 2nd Order ODE", xlabel="Time [days]", ylabel="Error [AU]")

# Add the errors for each method to the plot
for (method_label, (error, times, comp_time)) in results_dict
    plot!(p, times/day, error/AU, label="$(method_label) ($(comp_time) s)", linewidth=2)
end

# Display the plot
display(p)