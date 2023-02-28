using DifferentialEquations
using Plots
#blba

const GM = 1.327e20::Float64                            # G * M_sol [m^3 / s^2]
const AU = 1.496e11::Float64                            # 1AU [m]
const yr = 3.154e7::Float64                             # 1yr [s]
const day = 24.0 * 60.0 * 60.0::Float64                 # 1d [s]
const r_heliosphere = 100*AU                            # radius heliosphere [m]
v0 = 30e3                                               # speed of incomming isd in [m/s]
theta = pi/2
phi = 0

# parameters
tspan = (0.0,35*yr)                                    # time span (start, end)
dt = 1*day                                              # equidistant time step
no_of_trajectories = 50
s0  = [r_heliosphere * sin(theta)*cos(phi),r_heliosphere * sin(theta)*sin(phi), r_heliosphere * cos(theta),0,-v0,0]

#different entrypoints in heliosphere, acess via s0[:,i] for i 1:no_of_trajectories (0 kein gültiger index)]
for i in 1 : no_of_trajectories
    local phi =i / no_of_trajectories * pi
    if phi == pi/2
        continue                                        #skipp trajectory that would hit the sun
    end
    s = [r_heliosphere * sin(theta)*cos(phi),r_heliosphere * sin(theta)*sin(phi),r_heliosphere * cos(theta),0,-v0,0]
    global s0 = hcat(s0, s)
end

# Write the function (differential equation)
function EqOfMotion(ds, s, p, t, b=0)
    ds[1:3] = s[4:6]                                    # derivative of position = velocity
    ds[4:6] = grav(s[1:3]) - srp(s[1:3],b)              # derivative of velocity = acceleration
end

# Compute gravitational acceleration
function grav(r)
    dist = sqrt(r[1]^2 + r[2]^2 + r[3]^2)
    return - GM / dist^3 * r
end

# Compute solar radiation pressure acceleration
function srp(r,be)
    dist = sqrt(r[1]^2 + r[2]^2 + r[3]^2)
    return - be * GM / dist^3 * r
end

#solve for r
prob = [ODEProblem(EqOfMotion, s0[:,i], tspan) for i in 1:size(s0, 2)]                # ODE Problem
sol = [solve(prob[i], RK4(), adaptive=:false, dt = dt) for i in 1:size(prob, 1)]      # solve the problem
r = [Array{Float64}(undef, length(sol[i].t), 3) for i in 1:size(sol)[1]]          # prepare positions
for j in 1:size(r)[1]
    for i in eachindex(sol[j].t)                               # get positions
        r[j][i,1] = sol[j].u[i][1]
        r[j][i,2] = sol[j].u[i][2]
        r[j][i,3] = sol[j].u[i][3]
    end
end
  
r .= r ./ AU                                            # scale for plotting
size(r)[1]
#plot
plot(xlabel = "AU", ylabel = "AU")
scatter!([0],[0],label = "sun") # fügen das Sonnensymbol (Kreis mit Loch)
[plot!(r[i][:,1], r[i][:,2],label = false) for i in 1:size(r)[1]]
θ = LinRange(0 , 2*π , 100)
plot!(r_heliosphere./ AU  * cos.(θ),r_heliosphere./ AU  * sin.(θ),size=(400,400),label = "heliosphere")

