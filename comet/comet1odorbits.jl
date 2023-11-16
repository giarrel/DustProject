using CairoMakie
using LinearAlgebra

include("KeptoCartfkt.jl")
include("cometlaunch.jl")
include("comnames.jl")
include("constants.jl")

#method
method_string = [
    #RK/single step methods
    Euler(),
#    Midpoint(),
    ImplicitEuler(),
#    ImplicitMidpoint(),
#    Trapezoid(),
#    RK4(),
    Vern9(),


    #Adams/multistep methods
#    AB3(),
#    AB4(),
#    AB5(),
#    ABM32(),
#    ABM43(),
#    ABM54()
#

]


beta=0
stepsize=0.1day
periods = 1

comet_name = "Phaethon" 

# Initial conditions and problem setup
tspan = (0, periods * comets[comet_name].period)
initial_pos = keplerian_to_cartesian(comet_name, tspan[1], tspan[1])[1]
initial_vel = keplerian_to_cartesian(comet_name, tspan[1], tspan[1])[2]

# Initialize a dictionary to store errors, times, and computation times for each method
results_dict = Dict{String, Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}()

# Create a string to store the names of the used methods
methods_used = ""

# Iterate over the list of integration methods
for method in method_string
    
    local vel_num, pos_num, times = launch_from_comet_1ord(initial_vel, initial_pos, tspan, beta, method, stepsize=stepsize)

    local method_label = split(split(string(method), '{')[1], '(')[1]


    # Store the errors, times, and computation time for the current method
    results_dict[method_label] = (pos_num[1],pos_num[2],pos_num[3])

    # Add the name of the current method to the string
    global methods_used *= method_label * "_"
end

# Remove the trailing underscore from the methods_used string
methods_used = chop(methods_used)

f = Figure()
ax = Axis3(f[1, 1];title = "1st Order ODE Orbit step size = $(stepsize/day) days", xlabel="AU", ylabel="AU",zlabel="AU")

for (method_label, (pos_num)) in results_dict
    
    lines!(ax, pos_num[1]/AU,pos_num[2]/AU,pos_num[3]/AU, label="$(method_label)", linewith=1)
end

axislegend()

# Save the plot
save("DustProject/error_plots_newnew/1ode_$(comet_name)_stepsize$(stepsize/day)_days_periods$(periods)_methods_$(methods_used).pdf", f)
GC.gc()