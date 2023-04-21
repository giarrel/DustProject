#https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf

using Plots

include("KeptoCartfkt.jl")
include("cometlaunch.jl")

#constants
const GM = 1.327e20::Float64                            # G * M_sol [m^3 / s^2]
const AU = 1.496e11::Float64                            # 1AU [m]
const yr = 3.154e7::Float64                             # 1yr [s]
const day = 24.0 * 60.0 * 60.0::Float64                 # 1d [s]


# Orbital elements
e = 0.8898953806122173 # Eccentricity
a = 1.271341913087611 * AU # Semi-major axis (m)
i = 22.27353064493616 * (pi / 180) # Inclination (radians)
node = 265.1723386665869 * (pi / 180) # Longitude of ascending node (radians)
peri = 322.2166335445264 * (pi / 180) # Argument of perihelion (radians)
M0 = 196.5061141457209 * (pi / 180) # Mean anomaly (radians)
period = 523.5907604676752 * day # Period (s)
n = 0.6875598791667854 * (pi / 180) / day # Mean motion (radians/s)


# Generate orbit points
num_points = 1000
times = range(0, stop=period, length=num_points)
positions = [keplerian_to_cartesian(a, e, peri, node, i, M0, times[1], time)[1] for time in times]

# Extract X, Y, and Z components of positions
x_positions = [pos[1] for pos in positions]
y_positions = [pos[2] for pos in positions]
z_positions = [pos[3] for pos in positions]

# Plot orbit
plot3d(x_positions/AU, y_positions/AU, z_positions/AU, xlabel="X (AU)", ylabel="Y (AU)", zlabel="Z (AU)", title="Comet Orbit", label="from Kepler elements", aspect_ratio=1)

# Add sun at the center
scatter!([0], [0], [0], markersize=2, markercolor=:yellow, label="Sun")


initial_pos = keplerian_to_cartesian(a, e, peri, node, i, M0, times[1], times[1])[1]
initial_vel = keplerian_to_cartesian(a, e, peri, node, i, M0, times[1], times[1])[2]

# Initial conditions and problem setup
beta=0
tspan = (times[1], times[end])


vel_num ,pos_num = launch_from_comet(initial_vel, initial_pos, tspan, beta)

plot3d!(pos_num[1]/AU, pos_num[2]/AU, pos_num[3]/AU, label="from solver")
