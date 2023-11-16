using DifferentialEquations, GLMakie, LinearAlgebra
include("constants.jl")
include("duststartHS.jl")


Bfield_strength_after=1.5e-10
Bfield_strength_before=0.5e-10

function B_field(y_travel)
    y_travel_adjust=y_travel+80AU
    B=-(Bfield_strength_after-Bfield_strength_before)*tanh(100*y_travel_adjust)/2-(Bfield_strength_after-Bfield_strength_before)/2+Bfield_strength_after
    return B.*normalize([1,0,0]) 
    
end

function lorenza(q,m,v,y_travel)
    y_travel_v_jump=y_travel+80AU
    v_sol_jump = -(v_sol_cs-v_sol)*tanh(100*y_travel_v_jump)/2-(v_sol_cs-v_sol)/2+v_sol_cs
    return -q/m*cross((v-[0,v_sol_jump,0]),B_field(y_travel))
end


beta_params = [1]

r_dust_paramsm =[   1e-6
                    1e-7
                    1e-8
                    1e-9]

m_dust_paramsm=[    1e-14 1e-15 1e-16
                    1e-17 1e-18 1e-19
                    1e-20 1e-21 1e-22
                    1e-22 1e-23 1e-24]

turnparams = [  -0.3pi
                -0.2pi
                -0.1pi
                0
                0.1pi
                0.2pi
                0.3pi]

# Anfangsbedingungen
v00 = [0,26,0]* km/s # [vx(0), vy(0), vz(0)][0,26,0]* km/s 
x0 = [0, -81AU, 0.0] # [x(0), y(0), z(0)]


for i in eachindex(r_dust_paramsm)
    r_dust_params=r_dust_paramsm[i]
    m_dust_params=m_dust_paramsm[i,:]

for turn_param in turnparams

rot_mat = [1 0 0
            0 cos(turn_param) -sin(turn_param)
           0 sin(turn_param)  cos(turn_param)]
v0 = rot_mat*v00
startlocvels = [vcat(x0, v0)]

fig = Figure(resolution = (1.1*1920,1.1*1080))
#ax = Axis3(fig[1, 2], xlabel = "AU", ylabel = "AU", zlabel = "AU" ,title = "3d")
ax2 = Axis(fig[1, 1], xlabel = "AU", ylabel = "AU", title = "2d TS dust trajectories turned by $(turn_param/pi)pi")

globalmin=globalmax=0

for u0 in startlocvels
for m in m_dust_params
for r_dust in r_dust_params
for beta in beta_params

    local q=5*4*pi*eps_0*r_dust

    

    # Definition der Differentialgleichungen
    function charged_particle!(du, u, p, t)
        r = norm(u[1:3])
        du[1:3] = u[4:6]
        du[4:6] = lorenza(q,m,u[4:6],u[2])-GM * u[1:3] / r^3 * (1 - beta)
    end

    # Zeitbereich
    local tspan = (0, 200day)
    

    # Problem und Lösung
    local prob = ODEProblem(charged_particle!, u0, tspan)
    local sol = solve(prob,Vern9(),dt=0.001*day,adaptive=true,reltol=1e-20)

    localmax=maximum(Float32.(sol[3,:])/AU)
    localmin=minimum(Float32.(sol[3,:])/AU)
    if localmax > globalmax
        globalmax=localmax
    end
    if localmin < globalmin
        globalmin=localmin
    end



    # Positionen plotten
    #lines!(ax,label= string.(m), Float32.(sol[1,:])/AU, Float32.(sol[2,:])/AU, Float32.(sol[3,:])/AU, linewidth = 1)
    lines!(ax2,label= "mass = $(m)kg ,dust radius = $(r_dust)m , initial velocity = $(u0[4:6]/km)km/s , integrtaiontime=$((tspan[2]-tspan[1])/day)days , q/m=$(round(q/m,digits=3))c/kg",  Float32.(sol[2,:])/AU, Float32.(sol[3,:])/AU, linewidth = 1)
end
end
end
end

linex = LinRange(-80, -80, 10)
liney = LinRange(globalmin, globalmax, 10)


lines!(ax2, linex, liney, label="TS",linewidth = 3)

arrows!(ax2,[startlocvels[1][2]/AU],[startlocvels[1][3]/AU],[startlocvels[1][5]/100km],[startlocvels[1][6]/100km])



axislegend()
fig
save("C:\\Users\\lucac\\dustproject_clone\\DustProject\\TSplots\\TS_startloc_$(x0/AU)_size_$(r_dust_params)_turn$(turn_param/pi).png",fig)
end
end