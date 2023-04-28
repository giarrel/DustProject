using Plots
using LinearAlgebra

include("KeptoCartfkt.jl")
include("cometlaunch.jl")
include("comnames.jl")
include("constants.jl")


comet_name = "Phaethon"

# Initial conditions and problem setup
beta=0
tspan = (0, 5comets[comet_name].period)
method = Vern9()


# Add sun at the center
scatter([0], [0], [0], markersize=2, markercolor=:yellow, label="Sun")


initial_pos = keplerian_to_cartesian(comet_name, tspan[1], tspan[1])[1]
initial_vel = keplerian_to_cartesian(comet_name,tspan[1], tspan[1])[2]

vel_num ,pos_num , times= launch_from_comet(initial_vel, initial_pos, tspan, beta, method) 

#plot
method_label=split(split(string(method), '{')[1], '(')[1]
plot3d!(pos_num[1] ./ AU, pos_num[2] ./ AU, pos_num[3] ./ AU, label=method_label, aspect_ratio=1)#this can not handle to many values Arned one period and dt=100s doesnt work anymore

# Generate orbit points
positions = [keplerian_to_cartesian(comet_name, times[1], time)[1] for time in times]

# Extract X, Y, and Z components of positions
x_positions = [pos[1] for pos in positions]
y_positions = [pos[2] for pos in positions]
z_positions = [pos[3] for pos in positions]

# Plot orbit
plot3d!(x_positions/AU, y_positions/AU, z_positions/AU, xlabel="X (AU)", ylabel="Y (AU)", zlabel="Z (AU)", title="Comet Orbit", label="from Kepler elements")