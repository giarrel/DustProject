using DifferentialEquations
using LinearAlgebra
using Plots


const GM = 1.327e20::Float64                            # G * M_sol [m^3 / s^2]
const AU = 1.496e11::Float64                            # 1AU [m]
const yr = 3.154e7::Float64                             # 1yr [s]
const day = 24.0 * 60.0 * 60.0::Float64                 # 1d [s]
const r_heliosphere = 100*AU                            # radius heliosphere [m]
const SWSS = 2.5 *696340e3                              # 2.5*solar radius, [m]

# parameters
tspan = (0.0,60*yr)                                    # time span (start, end)
v0 = 30e3                                               # speed of incomming isd in [m/s]
beta=0.5                                                   #beta for srp
q_durch_m=0                                                #q/m for lorenz
integral_method ,relative_tollerance = RK4() , 1e-10           #in seconds
some_no_related_to_ammount_of_traj,dreide = 5,false     

function trajectories(ammount,dreid)
    s0 = [0,r_heliosphere,0,0,-v0,0]

    if dreid==false
        for i in 0:ammount
            s = [-r_heliosphere+2*r_heliosphere*i/ammount,r_heliosphere,0,0,-v0,0]
            s0=hcat(s0,s)
        end
        return s0
    else
        for i in 0 : ammount
            for j in 0 : ammount    
                s = [-r_heliosphere+2*r_heliosphere*i/ammount,r_heliosphere,-r_heliosphere+2*r_heliosphere*j/ammount,0,-v0,0]
                s0 = hcat(s0, s)
            end
        end
        return s0
    end
end

s0=trajectories(some_no_related_to_ammount_of_traj,dreide)

#use this for parameter studdy
#s0 = [AU,0,0,0,-v0,0]


# Write the function (differential equation)
function EqOfMotion(ds, s, p, t, b=beta,qm=q_durch_m)
    ds[1:3] = s[4:6]                                    # derivative of position = velocity
    ds[4:6] = grav_srp(s[1:3],b) +  lorenzf(qm,s[4:6],magnetic_field(s[1:3]))                    # derivative of velocity = acceleration
    #ds[4:6] = lorenzf(qm,s[4:6],magnetic_field(s[1:3])) #just to test what lorenzforce does
end

# Compute gravitational acceleration, Compute solar radiation pressure acceleration
function grav_srp(r,be)
    dist = sqrt(r[1]^2 + r[2]^2 + r[3]^2)
    return - GM / dist^3 * r + be * GM / dist^3 * r
end

function lorenzf(qm,v,B)
    return qm*cross(v,B)
end


# Definieren Sie die magnetische Feldfunktion der Parker Spirale
function magnetic_field(r_vec)
    # Konstanten
    B0 = 5e-9        # nT, interstellare Magnetfeldstärke bei 1 AU
    L = AU          # Längenskala, AU
    omega = 2.7e-6   # rad/s, Sonnenrotationsrate nach parkers paper oder 2π/(25.7*24*60*60)?


    # Umwandlung in zylindrische Koordinaten
    x, y, z = r_vec[1], r_vec[2], r_vec[3]
    r = sqrt(x^2 + y^2 + z^2)
    θ, φ =  acos(z/r), atan(y, x)
    v_r_solar = 400e3 #solar wind speed in m/s
    
    # Berechnen Sie die Parker-Spirale Magnetfeldkomponenten
    B_r = B0 * (L/r)^2
    B_phi = -B0 * omega * L^2 * sin(θ)/(v_r_solar*r)
    B_theta = 0.0
    
    # Umwandlung zurück in kartesische Koordinaten
    B_x, = B_r * sin(θ) * cos(φ) + B_theta * cos(θ) * cos(φ) - B_phi * sin(φ)
    B_y = B_r * sin(θ) * sin(φ) + B_theta * cos(θ) * sin(φ) + B_phi * cos(φ)
    B_z = B_r * cos(θ) - B_theta * sin(θ)

    return [B_x, B_y, B_z]
end

#callback stuff to end trajectory that hits the sun
function condition(dist,t,integrator) 
    s = integrator.u
    dist = sqrt(s[1]^2 + s[2]^2 + s[3]^2) # Die Entfernung vom Ursprung 
    dist < SWSS # Die Bedingung für das Stoppen end
end
function affect!(integrator)
    terminate!(integrator) # Die Aktion für das Stoppen
end

cb = DiscreteCallback(condition,affect!) # Erstellen Sie den Callback

#solve for r

prob = [ODEProblem(EqOfMotion, s0[:,i], tspan) for i in axes(s0,2)]                # ODE Problem
sol = [solve(prob[i], integral_method ,reltol = relative_tollerance, callback=cb) for i in axes(prob, 1)]      # solve the problem
r = [Array{Float64}(undef, length(sol[i].t), 3) for i in axes(sol)[1]]          # prepare positions
for j in axes(r)[1]
    for i in eachindex(sol[j].t)                               # get positions
        r[j][i,1] = sol[j].u[i][1]
        r[j][i,2] = sol[j].u[i][2]
        r[j][i,3] = sol[j].u[i][3]
    end
end

r .= r ./ AU                                            # scale for plotting

#plot2d function
function twodim()
    plot(xlabel = "AU", ylabel = "AU",xlims=(-110, 110), ylims=(-120, 140))
    scatter!([0],[0],label = "sun", color=:yellow) # fügen das Sonnensymbol (Kreis mit Loch)
    [plot!(r[i][:,1], r[i][:,2],label = false) for i in 1:size(r)[1]]
    θ = LinRange(0 , 2*π , 100)
    plot!(r_heliosphere./ AU  * cos.(θ),r_heliosphere./ AU  * sin.(θ),size=(400,400),label = "heliosphere",title="ISD Trajectories in 2D")
end

# 3D plot function
function threedim()
    plot3d()
    for i in eachindex(r)
        plot3d!(r[i][:,1], r[i][:,2], r[i][:,3],label = false)
    end

    print("last computed position of s0")
    print(r[1][size(r[1])[1],:])

    scatter3d!([0],[0],[0],label = "sun", color=:yellow)
    plot3d!(size=(800,600), xlabel="AU", ylabel="AU", zlabel="AU")
    plot3d!(title="ISD Trajectories in 3D")
end

#plots
#twodim()
threedim()