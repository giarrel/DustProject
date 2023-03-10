#looks bad

using Plots
using LinearAlgebra
using DifferentialEquations

#used Steyn (2020) https://iopscience.iop.org/article/10.3847/1538-4357/abb2a5/pdf

const AU = 1.496e11::Float64                            # 1AU [m]
const Rsun = 696340e3                                   # solar radius, [m]
const day = 24.0 * 60.0 * 60.0::Float64                 # 1d [s]

# Definieren Sie die magnetische Feldfunktion nach General-Fisk
function magnetic_field(r_vec)

    # Konstanten
    B0 = 5e-9        # nT, interstellare Magnetfeldstärke bei 1 AU
    #θ_ph_hm, θ_ss_mm, θ_ph_mm = 0.261799,pi/2,pi/2

    # Umwandlung in zylindrische Koordinaten
    x, y, z = r_vec[1], r_vec[2], r_vec[3]
    r = sqrt(x^2 + y^2 + z^2)
    θ, φ =  acos(z/r), atan(y, x)
    
    #benötigte werte für magnetfeld 
    OMEGA = 2.865e-6 #2*pi/(26*day)
    v_ref = 800e3                        #m/s
    v_solar = v_ref-OMEGA*r*sin(θ)   #solar wind speed in m/s from https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2002GL015028
    dv_dθ=-OMEGA*r*cos(θ)
    omega = fkt_omega(θ, φ)   # rad/s   the angular differential rotation rate
    phi0 = 0  #the azimuthal angle at t = 0
    phi_star = φ + OMEGA*r/v_solar-phi0
    #alpha = θ_ph_hm
    beta = 0.485#acos(1-(1-cos(θ_ss_mm))*(sin(alpha)^2/sin(θ_ph_mm)^2))-alpha

    
    # Magnetfeldkomponenten eq 12-14
    B_r = B0 * (1+(r*omega*sin(beta)*sin(phi_star))/(v_solar^2)*dv_dθ)

    B_theta = B_r*r/v_solar*omega*sin(beta)*sin(phi_star)

    B_phi = B_r*r/v_solar*  (omega*sin(beta)*cos(θ)*cos(phi_star)+
                            omega*cos(beta)*sin(θ)-OMEGA*sin(θ))
    
    
    # Umwandlung zurück in kartesische Koordinaten
    B_x, = B_r * sin(θ) * cos(φ) + B_theta * cos(θ) * cos(φ) - B_phi * sin(φ)
    B_y = B_r * sin(θ) * sin(φ) + B_theta * cos(θ) * sin(φ) + B_phi * cos(φ)
    B_z = B_r * cos(θ) - B_theta * sin(θ)

    return [B_x, B_y, B_z]
end

function fkt_omega(θ, φ)
    omega = 0.464 * 1e-6 * cos(θ)^2 + 0.328 * 1e-6 *cos(θ)^4
    return omega
end

#magn field lines calculated and plotted from field

plot3d(size=(800,800,800))
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

timespan=(0.0,1e14)
θ , φ = pi/2 , 0                            #θ 0 bis pi, φ 0 bis 2pi
startpoint= 2.5*Rsun .* [sin(θ)*cos(φ),sin(θ)*sin(φ),cos(θ)]

magn_field_line(startpoint,timespan,θ,φ)