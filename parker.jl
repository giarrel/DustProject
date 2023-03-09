using Plots
#using CairoMakie
using LinearAlgebra
using DifferentialEquations

#used Kinematic models of the interplanetary magnetic field Christoph Lhotka and Yasuhito Narita (2019), The Heliospheric Magnetic Field Mathew J. Owens(2013), nicht in datenbank gefunden hat aber geholfen die beiden anderen paper zu finden : Heliospheric Magnetic Field and The Parker Model N. S. Svirzhevsky(2021) 

const AU = 1.496e11::Float64                            # 1AU [m]

# Definieren Sie die magnetische Feldfunktion der Parker Spirale
function magnetic_field(r_vec, t = 0)
    # Konstanten
    B0 = 5e-9        # nT, interstellare Magnetfeldstärke bei 1 AU
    L = AU          # Längenskala, AU
    omega = 2.7e-6   # rad/s, Sonnenrotationsrate nach parkers paper oder 2π/(25.7*24*60*60)?


    # Umwandlung in zylindrische Koordinaten
    x, y, z = r_vec[1], r_vec[2], r_vec[3]
    r = sqrt(x^2 + y^2 + z^2)
    θ, φ =  acos(z/r), atan(y, x)+ omega*t
    v_r_solar = 400e3 #solar wind speed in m/s
    
    # Berechnen Sie die Parker-Spirale Magnetfeldkomponenten neu: cos(teta) eingefügt nach parkers paper rauskomentiert wegen nicht sicher
    B_r = B0 * (r/L)^(-2) #* cos(θ)
    B_phi = -B0 * omega * L^2 * sin(θ)/(v_r_solar*r) #* cos(θ)
    B_theta = 0.0
    
    # Umwandlung zurück in kartesische Koordinaten
    B_x, = B_r * sin(θ) * cos(φ) + B_theta * cos(θ) * cos(φ) - B_phi * sin(φ)
    B_y = B_r * sin(θ) * sin(φ) + B_theta * cos(θ) * sin(φ) + B_phi * cos(φ)
    B_z = B_r * cos(θ) - B_theta * sin(θ)

    return [B_x, B_y, B_z]
end

#magn field lines calculated and plotted from field
plot3d(xlims=(-100,100),ylims=(-100,100),zlims=(-100,100),size=(800,800,800))

function ode_system(du,u,p,t)
    du .= magnetic_field(u) ./ norm(magnetic_field(u))
end

function magn_field_line(u0, tspan) #need to solve g' = F(g)
    prob = ODEProblem(ode_system,u0,tspan)
    sol = solve(prob,dtmax=1e10)
    X = [x[1] for x in sol.u]
    Y = [x[2] for x in sol.u]
    Z = [x[3] for x in sol.u]
    u0AU=u0./AU
    plot3d!(X/AU,Y/AU,Z/AU,label="Parker Spiral with starting at $u0AU AU")
    xlabel!("AU")
    ylabel!("AU")
    zlabel!("AU")
    
end

timespan=(0.0,1e15)
startpoint=[AU/10,AU/10,AU/10]

magn_field_line(startpoint,timespan)

#=
#first approach using quiver maybee need later to see how i created a grid

xmin, xmax = -10*AU, 10*AU
ymin, ymax = -10*AU, 10*AU

# Erstellen eines Gitters von Koordinaten im Raum, mit entsprechenen magnetfeld werten
xgrid, ygrid = range(xmin, xmax, length=20), range(ymin, ymax, length=20)

grid=[[x,y] for x in xgrid for y in ygrid] #vektor länge lenth *length mit einträgen [x,y] aus allen x y kombinationen möglichkeiten
X,Y = getindex.(grid, 1),getindex.(grid, 2) #erste einträge von vektor in vektor, zweite einträge von vektor in vektor

B = [1e9*magnetic_field([x,y,0],1/2)[1:2] for x in xgrid for y in ygrid]#analog zu grid
Bx,By=getindex.(B, 1),getindex.(B, 2)#analog zu X,Y


#f(x,y)=Point2f(magnetic_field([x,y,0])[1:2])
#norm(f(10*AU,0))

plot(size=(400,400),xlims=(-10,10),ylims=(-10,10),title = "Parker Spiral z=0")
quiver!(X/AU,Y/AU,quiver=(Bx,By)) #plot(X/AU,Y/AU) mit pfeilen (Bx,By)  (geht nur 2d)
xlabel!("AU")
ylabel!("AU")
#p=streamplot(f,xmin..xmax,ymin..ymax)
#display(p)
=#