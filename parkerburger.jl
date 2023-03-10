#ACHTUNG NOCH nicht fertig/lauffähig

using Plots
using LinearAlgebra
using DifferentialEquations

#used Burger (2008) https://iopscience.iop.org/article/10.1086/525039/pdf

const AU = 1.496e11::Float64                            # 1AU [m]
const Rsun = 696340e3                                   # solar radius, [m]

# Definieren Sie die magnetische Feldfunktion der Parker Spirale
function magnetic_field(r_vec, t = 0)
    # Konstanten
    B0 = 5e-9        # nT, interstellare Magnetfeldstärke bei 1 AU
    L = AU          # Längenskala, AU
    omega = 2.865e-6   # rad/s, Sonnenrotationsrate nach parkers paper oder 2π/(25.7*24*60*60)?


    # Umwandlung in zylindrische Koordinaten
    x, y, z = r_vec[1], r_vec[2], r_vec[3]
    r = sqrt(x^2 + y^2 + z^2)
    θ, φ =  acos(z/r), atan(y, x)+ omega*t
    v_r_solar = 430e3 #solar wind speed in m/s
    
    # Berechnen Sie die Parker-Spirale Magnetfeldkomponenten von Burger eq 4
    B_r = B0 * (L/r)^2

    B_theta = B_r*r/v_r_solar*omega_star*sin(beta_star)*sin(phi_star)

    B_phi = B_r*r/v_r_solar*(omega_star*sin(beta_star)*cos(θ)*cos(phi_star)+
                            sin(θ)*(omega_star*cos(beta_star)-omega)+
                            domega_star_dtheta*sin(beta_star)*sin(θ)*cos(phi_star)+
                            omega_star*dbeta_star_dtheta*cos(beta_star)*sin(θ)*cos(phi_star))
    
    
    # Umwandlung zurück in kartesische Koordinaten
    B_x, = B_r * sin(θ) * cos(φ) + B_theta * cos(θ) * cos(φ) - B_phi * sin(φ)
    B_y = B_r * sin(θ) * sin(φ) + B_theta * cos(θ) * sin(φ) + B_phi * cos(φ)
    B_z = B_r * cos(θ) - B_theta * sin(θ)

    return [B_x, B_y, B_z]
end

#magn field lines calculated and plotted from field

plot3d(xlims=(-100,100),ylims=(-100,100),zlims=(-100,100),size=(800,800,800))
scatter3d!([0,0,0],[0,0,0],[0,0,0],label="sun", color=:yellow)

function ode_system(du,u,p,t)
    du .= magnetic_field(u) ./ norm(magnetic_field(u))
end

function magn_field_line(u0, tspan, θ , φ) #need to solve g' = F(g) siehe https://uwaterloo.ca/physics-of-information-lab/sites/ca.physics-of-information-lab/files/uploads/files/amath231_1.2.2.pdf
    prob = ODEProblem(ode_system,u0,tspan)
    sol = solve(prob,dtmax=1e10)
    X = [x[1] for x in sol.u]
    Y = [x[2] for x in sol.u]
    Z = [x[3] for x in sol.u]
    θpi=θ/π
    φpi=φ/π
    plot3d!(X/AU,Y/AU,Z/AU,label="Parker Spiral with starting at θ=$θpi π φ=$φpi π")
    xlabel!("AU")
    ylabel!("AU")
    zlabel!("AU")
    
end

timespan=(0.0,1e15)
θ , φ = π/2 , 0                            #θ 0 bis pi, φ 0 bis 2pi
startpoint= Rsun .* [sin(θ)*cos(φ),sin(θ)*sin(φ),cos(θ)]

magn_field_line(startpoint,timespan,θ,φ)