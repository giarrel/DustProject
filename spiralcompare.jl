using Plots
using LinearAlgebra
using DifferentialEquations

#used Kinematic models of the interplanetary magnetic field Christoph Lhotka and Yasuhito Narita (2019), The Heliospheric Magnetic Field Mathew J. Owens(2013), nicht in datenbank gefunden hat aber geholfen die beiden anderen paper zu finden : Heliospheric Magnetic Field and The Parker Model N. S. Svirzhevsky(2021) 

const AU = 1.496e11::Float64                            # 1AU [m]
const Rsun = 696340e3                                   # solar radius, [m]

# Definieren Sie die magnetische Feldfunktion der Parker Spirale
function magnetic_field(r_vec)
    # Konstanten
    B0 = 5e-9        # nT, interstellare Magnetfeldstärke bei 1 AU
    L = AU          # Längenskala, AU
    omega = 2.865e-6   # rad/s, Sonnenrotationsrate nach parkers paper oder 2π/(25.7*24*60*60)?


    # Umwandlung in zylindrische Koordinaten
    x, y, z = r_vec[1], r_vec[2], r_vec[3]
    r = sqrt(x^2 + y^2 + z^2)
    θ, φ =  acos(z/r), atan(y, x)
    v_r_solar = 430e3 #solar wind speed in m/s
    
    # Berechnen Sie die Parker-Spirale Magnetfeldkomponenten neu: cos(teta) eingefügt nach parkers paper rauskomentiert wegen nicht sicher, sollte passen siehe kleinmann bei eq 7
    B_r = B0 * (r/L)^(-2) * cos(θ)/abs(cos(θ))
    B_phi = -B0 * omega * L^2 * sin(θ)/(v_r_solar*r) * cos(θ)/abs(cos(θ))
    B_theta = 0.0
    
    # Umwandlung zurück in kartesische Koordinaten
    B_x, = B_r * sin(θ) * cos(φ) + B_theta * cos(θ) * cos(φ) - B_phi * sin(φ)
    B_y = B_r * sin(θ) * sin(φ) + B_theta * cos(θ) * sin(φ) + B_phi * cos(φ)
    B_z = B_r * cos(θ) - B_theta * sin(θ)

    return [B_x, B_y, B_z]
end

# Funktion zum Hinzufügen der analytischen Parker-Spirale zum 3D-Plot
function add_parker_spiral_to_plot!(plt)
    rmin = 1.0 * Rsun  # inner radius, [m]
    rmax = 5.3 * AU    # outer radius, [m]

    r, φ = parker_spiral_2d(rmin, rmax)
    X = r .* cos.(φ) ./ AU
    Y = r .* sin.(φ) ./ AU
    Z = zeros(length(r))

    plot!(plt, X, Y, Z, label="Analytical Parker Spiral", color=:green, lw=2)
end

# Erstellen Sie den 3D-Plot
#magn field lines calculated and plotted from field

plot3d(xlims=(-1,1),ylims=(-1,1),size=(800,800,800))
scatter3d!([0,0,0],[0,0,0],[0,0,0],label="sun", color=:yellow)

function ode_system(du,u,p,t)
    du .= magnetic_field(u) ./ norm(magnetic_field(u))
end

function magn_field_line(u0, tspan, θ , φ) #need to solve g' = F(g) siehe https://uwaterloo.ca/physics-of-information-lab/sites/ca.physics-of-information-lab/files/uploads/files/amath231_1.2.2.pdf
    prob = ODEProblem(ode_system,u0,tspan)
    sol = solve(prob,reltol=1e-17)
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
startpoint= 2.5*Rsun .* [sin(θ)*cos(φ),sin(θ)*sin(φ),cos(θ)]

magn_field_line(startpoint,timespan,θ,φ)

# Fügen Sie die analytische Parker-Spirale hinzu
add_parker_spiral_to_plot!(plt)

# Zeigen Sie den Plot an
display(plt)