using DifferentialEquations
using Plots


const GM = 1.327e20::Float64                            # G * M_sol [m^3 / s^2]
const AU = 1.496e11::Float64                            # 1AU [m]
const yr = 3.154e7::Float64                             # 1yr [s]
const day = 24.0 * 60.0 * 60.0::Float64                 # 1d [s]
const r_heliosphere = 100*AU                            # radius heliosphere [m]
v0 = 26e3                                               # speed of incomming isd in [m/s]
theta = pi/2
phi = 0

# parameters
tspan = (0.0,40*yr)                                    # time span (start, end)
dt = 1*day                                              # equidistant time step
no_of_trajectories = 10
#s0  = [r_heliosphere * sin(theta)*cos(phi),r_heliosphere * sin(theta)*sin(phi), r_heliosphere * cos(theta),0,-v0,0]
s0 = [-r_heliosphere+2*r_heliosphere,r_heliosphere,0,0,-v0,0]

#different entrypoints in heliosphere, acess via s0[:,i] for i 1:no_of_trajectories (0 kein gültiger index)]
#=for i in 1 : no_of_trajectories
    local phi =i / no_of_trajectories * pi
    if phi == pi/2
        continue                                        #skipp trajectory that would hit the sun
    end
    s = [r_heliosphere * sin(theta)*cos(phi),r_heliosphere * sin(theta)*sin(phi),r_heliosphere * cos(theta),0,-v0,0]
    global s0 = hcat(s0, s)
end=#

for i in 1 : no_of_trajectories
    for j in 0 : no_of_trajectories    
        s = [-r_heliosphere+2*r_heliosphere*i/no_of_trajectories,r_heliosphere,-r_heliosphere+2*r_heliosphere*j/no_of_trajectories,0,-v0,0]
        global s0 = hcat(s0, s)
    end
end


# Write the function (differential equation)
function EqOfMotion(ds, s, p, t, b=10)
    ds[1:3] = s[4:6]                                    # derivative of position = velocity
    ds[4:6] = grav_srp(s[1:3],b)                        # derivative of velocity = acceleration
end

# Compute gravitational acceleration, Compute solar radiation pressure acceleration
function grav_srp(r,be)
    dist = sqrt(r[1]^2 + r[2]^2 + r[3]^2)
    return - GM / dist^3 * r + be * GM / dist^3 * r
end


function condition(dist,t,integrator) 
    s = integrator.u
    dist = sqrt(s[1]^2 + s[2]^2 + s[3]^2) # Die Entfernung vom Ursprung 
    dist <0.1*AU # Die Bedingung für das Stoppen end
end
function affect!(integrator)
    terminate!(integrator) # Die Aktion für das Stoppen
end

cb = DiscreteCallback(condition,affect!) # Erstellen Sie den Callback


#solve for r
prob = [ODEProblem(EqOfMotion, s0[:,i], tspan) for i in 1:size(s0, 2)]                # ODE Problem
sol = [solve(prob[i], RK4(), adaptive=:false, dt = dt,callback=cb) for i in 1:size(prob, 1)]      # solve the problem
r = [Array{Float64}(undef, length(sol[i].t), 3) for i in 1:size(sol)[1]]          # prepare positions
for j in 1:size(r)[1]
    for i in eachindex(sol[j].t)                               # get positions
        r[j][i,1] = sol[j].u[i][1]
        r[j][i,2] = sol[j].u[i][2]
        r[j][i,3] = sol[j].u[i][3]
    end
end
  
r .= r ./ AU                                            # scale for plotting

#plot
plot(xlabel = "AU", ylabel = "AU",xlims=(-110, 110), ylims=(-120, 140))
scatter!([0],[0],label = "sun") # fügen das Sonnensymbol (Kreis mit Loch)
[plot!(r[i][:,1], r[i][:,2],label = false) for i in 1:size(r)[1]]
θ = LinRange(0 , 2*π , 100)
plot!(r_heliosphere./ AU  * cos.(θ),r_heliosphere./ AU  * sin.(θ),size=(400,400),label = "heliosphere",title="ISD Trajectories in 2D")
#
# 3D plot
plot3d(r[1][:,1], r[1][:,2], r[1][:,3],label = false)
for i in 2:length(r)
    plot3d!(r[i][:,1], r[i][:,2], r[i][:,3],label = false)
end

scatter3d!([0],[0],[0],label = "sun", markersize=10)

plot3d!(size=(800,600), xlabel="AU", ylabel="AU", zlabel="AU")
plot3d!(title="ISD Trajectories in 3D")
#
