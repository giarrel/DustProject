using DifferentialEquations, GLMakie, LinearAlgebra
include("constants.jl")
include("duststartHS.jl")


d = 1/2*sol_day*v_sol_cs #distance between 2 current sheetsurfaces
Bfield_strength=1.5e-10


function B_field(y_travel,t)

    y_travel_rel=y_travel-t*v_sol_cs #Field and CS moving with Solar wind

    if iseven(div(y_travel_rel,d)) == true #bestimmt vorzeichen von magnetfeld
        B=Bfield_strength
    else
        B=-Bfield_strength
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

turnparams = [
                
                0
                
                ]

tstartvalsv = [ 
    
    d/v_sol_cs-d/v_sol_cs*3/2
]


# Anfangsbedingungen
v00 = [0,26,0]* km/s # [vx(0), vy(0), vz(0)][0,26,0]* km/s
x0 = [0, 130d, 0.0] # [x(0), y(0), z(0)]

for i in eachindex(r_dust_paramsm)
    r_dust_params=r_dust_paramsm[i]
    m_dust_params=m_dust_paramsm[i,:]

for turn_param in turnparams

rot_mat = [1 0 0
            0 cos(turn_param) -sin(turn_param)
           0 sin(turn_param)  cos(turn_param)]
v0 = rot_mat*v00
startlocvels = [vcat(x0, v0)]

for tstartvals in tstartvalsv


fig = Figure(resolution = (1.1*1920,1.1*1080))
#ax = Axis3(fig[1, 1], xlabel = "AU", ylabel = "AU", zlabel = "AU" ,title = "3d")
ax2 = Axis(fig[1, 1], xlabel = "AU", ylabel = "AU", title = "2d HS solar maximum trajectories turned by $(turn_param/pi)pi with stepsize 0.0005day")

for u0 in startlocvels
for tstart in tstartvals
for m in m_dust_params
for r_dust in r_dust_params
for beta in beta_params

    
    local q=5*4*pi*eps_0*r_dust


    function lorenza(q,m,v,y_travel,t)
        return -q/m*cross((v-[0,v_sol_cs,0]),B_field(y_travel,t))#-[0,v_sol_cs,0]
    end

    

    # Definition der Differentialgleichungen
    function charged_particle!(du, u, p, t)
        r = norm(u[1:3])
        du[1:3] = u[4:6]
        du[4:6] = lorenza(q,m,u[4:6],u[2],t)#-GM * u[1:3] / r^3 * (1 - beta)
    end

    # Zeitbereich
    local tspan = (tstart, 600day)
    

    # Problem und Lösung
    local prob = ODEProblem(charged_particle!, u0, tspan)
    local sol = solve(prob,Vern9(),dt=0.0005*day,adaptive=false,reltol=1e-20)

    t1=t2=told=0

    for i in 1:(length(sol.t)-1)
        

        if B_field(sol.u[i][2],sol.t[i])[1]>0
            if B_field(sol.u[i+1][2],sol.t[i+1])[1]<0
                t1=sol.t[i]
                #print(sol.t[i])
                #print("---")
            end
        end


        if B_field(sol.u[i][2],sol.t[i])[1]<0
            if B_field(sol.u[i+1][2],sol.t[i+1])[1]>0
                t2=sol.t[i]
                #print(sol.t[i])
                #print("---")
            end
        end


        if (t1-t2) == told
            nothing
        else
            #print((t1-t2),"^^^^")#prints time between polarities
        end

        told=(t1-t2)
        #print(norm(sol.u[i][4:6]),"---")
    end


    # Positionen plotten
    #lines!(ax,label= string.(m), Float32.(sol[1,:])/AU, Float32.(sol[2,:])/AU, Float32.(sol[3,:])/AU, linewidth = 1)
    lines!(ax2,label= "m=$(m)kg ,r=$(r_dust)m ,v_y=$(u0[5]/km)km/s ,v_z=$(u0[6]/km)km/s ,start delay= $(tstart/day)days, integrtaiontime=$(sol.t[end]/day)days, q/m=$(round(q/m,digits=3))c/kg",  Float32.(sol[2,:])/AU, Float32.(sol[3,:])/AU, linewidth = 1)
    
end
end
end
end
end

axislegend()
fig
save("C:\\Users\\lucac\\dustproject_clone\\DustProject\\HSintegratorplots\\HS_time_$(tstartvals/d*v_sol_cs)_size_$(r_dust_params)_turn$(turn_param/pi)pi_stepsize_0.0005.png",fig)
end
end
end