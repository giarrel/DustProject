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
tspan = (0.0,2yr)                                    # time span (start, end)
v0 = 5e3                                               # speed of incomming isd in [m/s]
beta=0                                                   #beta for srp
q_durch_m=0                                                #q/m for lorenz
fixed_step_size = day                                    #in seconds  
s0 = [AU,0,0,0,-v0,0]

# Definieren Sie verschiedene Integrationsmethoden und Farben für den Plot
integration_methods = [Tsit5(), BS3(), Rodas5P()]#, DP5(), SSPRK22(),Rodas5P()]
method_names = [string(typeof(method)) for method in integration_methods]
method_labels = [split(string(method), '(')[1] for method in integration_methods]


# Write the function (differential equation)
function EqOfMotion(ds, s, p, t, b=beta,qm=q_durch_m)
    ds[1:3] = s[4:6]                                    # derivative of position = velocity
    ds[4:6] = grav_srp(s[1:3],b) +  lorenzf(qm,s[4:6],magnetic_field(s[1:3]))                    # derivative of velocity = acceleration
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

# Nach der Vereinfachung:
prob = ODEProblem(EqOfMotion, s0, tspan)
solutions = [solve(prob, method, dt=fixed_step_size, callback=cb) for method in integration_methods]

# In der 3D-Plot-Funktion:
function threedim()
    plot3d()
    for (method_index, sol) in enumerate(solutions)
        r = Array{Float64}(undef, length(sol.t), 3)
        for j in eachindex(sol.t)
            r[j, 1] = sol.u[j][1]
            r[j, 2] = sol.u[j][2]
            r[j, 3] = sol.u[j][3]
        end
        r .= r ./ AU
        plot3d!(r[:, 1], r[:, 2], r[:, 3], label = method_labels[method_index])
    end

    scatter3d!([0], [0], [0], label = "sun", color = :yellow)
    plot3d!(size = (800, 600), xlabel = "AU", ylabel = "AU", zlabel = "AU")
    plot3d!(title = "Test Trajectories in 3D with dt=$fixed_step_size seconds and v0=$v0 m/s")
    plot3d!(legend = :outertopright)
end

# Plots
threedim()