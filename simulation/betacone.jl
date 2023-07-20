using DifferentialEquations
using LinearAlgebra
using Plots
using PlotlyJS

include("constants.jl")
include("duststartforbeta.jl")

#paramerters
tspan=20yr
beta=1.05
vd = 26.0e3 # ISD dust speed at infinity [SI]
no_of_starts=1
u0=startparams(no_of_starts,vd,beta)
p = (beta, GM)
au=AU

function trajectory_ode(du, u, p, t)
        
    beta, GM = p
    r = norm(u[1:3])
    du[1:3] = u[4:6]
    du[4:6] = -GM * u[1:3] / r^3 * (1 - beta)

end

#draw beta exclusion zones

x_pol(r,phi) = r>=0 ? r * cosd(phi) : 0 #convert polar coordinates to Cartesian

y_pol(r,phi) = r>=0 ? r * sind(phi) : 0 #convert polar coordinates to Cartesian

R_excl1(beta, phi) = -4*GM*(1-beta)/vd^2 ./ (1.0.+cosd.(phi)) # beta exclusion parabola
#derivative(beta, phi) = -4*GM*(1-beta)/vd^2 ./ ((1.0.+cosd.(phi)).*(1.0.+cosd.(phi))) .*sind.(phi)
#tangente(phi) = R_excl1(beta, phi) + derivative(beta, phi)*()

exclusionzone(beta, phi)=(x_pol.(R_excl1(beta, phi),phi), y_pol.(R_excl1(beta, phi),phi))

phi = -179:179 # all angles

R1 = R_excl1(beta, phi)/au
# parabola radius in dependence of angle
plotlyjs()
plo=plot(x_pol.(R1,phi), y_pol.(R1,phi), xlims=(-10,10), ylims=(-5,5),label="β = $beta",size=(1200,800))


for i in axes(u0)[2]
    local prob = ODEProblem(trajectory_ode, u0[1:6,i], tspan, p)
    local sol = solve(prob, Vern9(),reltol=1e-15)#  adaptive=false, dt=0.01day)#

    # Extract the results
    local x_pos_numeric = [point[1] for point in sol.u]
    local y_pos_numeric = [point[2] for point in sol.u]
    local z_pos_numeric = [point[3] for point in sol.u]

    min_distances=zeros(length(x_pos_numeric))

    for j in axes(x_pos_numeric)[1]
        min_dist(phi) = norm([x_pos_numeric[j]-exclusionzone(beta,phi)[1],y_pos_numeric[j]-exclusionzone(beta,phi)[2]])
        min_distances[j],_ = findmin(min_dist,-179:179)
    end

    min_val,min_index=findmin(min_distances)
    
    alpha=acosd(x_pos_numeric[min_index]/norm([x_pos_numeric[min_index],y_pos_numeric[min_index]]))

    println("relativer fehler min abstand zum cone: ",(norm([x_pos_numeric[min_index],y_pos_numeric[min_index]])-R_excl1(beta, alpha))/R_excl1(beta, alpha))

    plot!(plo,x_pos_numeric/AU, y_pos_numeric/AU)
end

display(plo)
GC.gc()