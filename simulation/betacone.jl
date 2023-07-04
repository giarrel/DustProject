using DifferentialEquations
using LinearAlgebra
using Plots

include("constants.jl")
include("duststart.jl")

#paramerters
tspan=20yr
beta=1.5
u0=startparams(10,26.0e3/(1e3),0,0)
p = (beta, GM)
au=AU

function trajectory_ode(du, u, p, t)
        
    beta, GM = p
    r = norm(u[1:3])
    du[1:3] = u[4:6]
    du[4:6] = -GM * u[1:3] / r^3 * (1 - beta)

end

#draw beta exclusion zones

vd = 26.0e3 # ISD dust speed at infinity [SI]

x_pol(r,phi) = r>=0 ?
r * cosd(phi) : 0 #convert polar coordinates to Cartesian

y_pol(r,phi) = r>=0 ?
r * sind(phi) : 0 #convert polar coordinates to Cartesian

R_excl1(beta, phi) = -4*GM*(1-beta)/vd^2 ./ (1.0.+cosd.(phi)) # beta exclusion parabola

phi = 0:360 # all angles

R1 = R_excl1(1.5, phi)/au
# parabola radius in dependence of angle

plo=plot(x_pol.(R1,phi), y_pol.(R1,phi), xlims=(-5.5,5.5), ylims=(-5.5,5.5),label="β = 1.5")

for i in axes(u0)[2]
    local prob = ODEProblem(trajectory_ode, u0[1:6,i], tspan, p)
    local sol = solve(prob, RadauIIA5(), adaptive=false, dt=0.01day)

    # Extract the results
    local x_pos_numeric = [point[1] for point in sol.u]
    local y_pos_numeric = [point[2] for point in sol.u]
    local z_pos_numeric = [point[3] for point in sol.u]
    
    plot!(plo,x_pos_numeric/AU, y_pos_numeric/AU)
end

display(plo)
GC.gc()