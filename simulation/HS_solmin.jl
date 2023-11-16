using DifferentialEquations, GLMakie, LinearAlgebra
include("constants.jl")
include("duststartHS.jl")

d = 1/2*sol_day*v_sol_cs #distance between 2 current sheetsurfaces
Bfield_strength=1.5e-10

Bvorzeichen = [1
-1]


function B_field_min(z)

    #B=-2*Bfield_strength*tanh(1000*z)/2  # no much difference model breakks anyways when too close to zero

    
    if z>0
        B=-Bfield_strength
    else
        B=Bfield_strength
    end
    

    return [B,0,0]
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

turnparams = [  -0.1pi
-0.05pi
0
0.05pi
0.1pi]

tstartvalsv = [ 
            0
                ]

# Anfangsbedingungen
v00 = [0,26,0]* km/s # [vx(0), vy(0), vz(0)][0,26,0]* km/s
x0 = [0, 130d, 0.1*AU] # [x(0), y(0), z(0)]


for i in eachindex(r_dust_paramsm)
    r_dust_params=r_dust_paramsm[i]
    m_dust_params=m_dust_paramsm[i,:]

for turn_param in turnparams

rot_mat = [ 1 0 0
            0 cos(turn_param) -sin(turn_param)
            0 sin(turn_param)  cos(turn_param)]
v0 = rot_mat*v00
startlocvels = [vcat(x0, v0)]

for tstartvals in tstartvalsv
for vorzeichen in Bvorzeichen

fig = Figure(resolution = (1.1*1920,1.1*1080))
#ax = Axis3(fig[1, 1], xlabel = "AU", ylabel = "AU", zlabel = "AU" ,title = "3d at sol min")
ax2 = Axis(fig[1, 1], xlabel = "AU", ylabel = "AU", title = "2d HS trajectories at at sol min")

globalmin=globalmax=-90

for u0 in startlocvels
for tstart in tstartvals
for m in m_dust_params
for r_dust in r_dust_params
for beta in beta_params

    local q=5*4*pi*eps_0*r_dust

    function lorenza(q,m,v,z_travel)
        return -q/m*cross((v-[0,v_sol_cs,0]),vorzeichen*B_field_min(z_travel))#-[0,v_sol_cs,0]
    end

    

    # Definition der Differentialgleichungen
    function charged_particle!(du, u, p, t)
        r = norm(u[1:3])
        du[1:3] = u[4:6]
        du[4:6] = lorenza(q,m,u[4:6],u[3])-GM * u[1:3] / r^3 * (1 - beta)
    end

    # Zeitbereich
    local tspan = (tstart, 1200day)
    

    # Problem und Lösung
    local prob = ODEProblem(charged_particle!, u0, tspan)
    local sol = solve(prob,Vern9(),dt=0.001*day,adaptive=true,reltol=1e-20)
    t1=t2=told=0

    for i in 1:(length(sol.t)-1)
        nothing #can check stuff in sols
    end

    localmax=maximum(Float32.(sol[2,:])/AU)
    localmin=minimum(Float32.(sol[2,:])/AU)
    if localmax > globalmax
        globalmax=localmax
    end
    if localmin < globalmin
        globalmin=localmin
    end

    # Positionen plotten
    #lines!(ax,label= string.(m), Float32.(sol[1,:])/AU, Float32.(sol[2,:])/AU, Float32.(sol[3,:])/AU, linewidth = 1)
    lines!(ax2,label= "m=$(m)kg ,r=$(r_dust)m ,v_y=$(round(u0[5]/km,digits=2))km/s ,v_z=$(round(u0[6]/km,digits=2))km/s , integrtaiontime=$(sol.t[end]/day)days",  Float32.(sol[2,:])/AU, Float32.(sol[3,:])/AU, linewidth = 1)
    
end
end
end
end
end
liney = LinRange(0, 0, 10)
linex = LinRange(globalmin, globalmax, 10)

lines!(ax2, linex, liney, label="CS",linewidth = 3)
arrows!(ax2,[startlocvels[1][2]/AU],[startlocvels[1][3]/AU],[startlocvels[1][5]/100km],[startlocvels[1][6]/100km],normalize=true)
axislegend()
fig
save("C:\\Users\\lucac\\dustproject_clone\\DustProject\\HSminplots\\HSmin_time_$(tstartvals/d*v_sol_cs)_size_$(r_dust_params)_turn$(turn_param/pi)_vorzeichen$(vorzeichen).png",fig)
end
end
end
end