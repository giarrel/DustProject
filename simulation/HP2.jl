using DifferentialEquations, GLMakie, LinearAlgebra
include("constants.jl")
include("duststartHS.jl")


Bfield_strength_after=6.0e-10
Bfield_strength_before=1.5e-10
turndistoutsidem = [0
                    0.5AU
                    1AU
                    2AU
                    3AU
                    5AU
                    1e10AU]

turnparams = [
                
                -0.25pi
                -0.12pi
                -0.05pi
                0
                0.05pi
                0.12pi
                0.25pi
                
]

random_turn_param=0.0*2pi
rot_matB = [ cos(random_turn_param) 0 -sin(random_turn_param)
            0                      1        0
            sin(random_turn_param) 0 cos(random_turn_param)]

B_direction_inside=[1,0,0]
B_direction_outside=*(rot_matB,B_direction_inside)


function B_field(y_travel)
    y_travel_adjust=y_travel+100AU
    B=-(Bfield_strength_after-Bfield_strength_before)*tanh(100*y_travel_adjust)/2-(Bfield_strength_after-Bfield_strength_before)/2+Bfield_strength_after

    if y_travel>-100AU
        return B.*B_direction_inside
    else
        return B.*B_direction_outside
    end
end


# Anfangsbedingungen
v00 = [0,26,0]* km/s # [vx(0), vy(0), vz(0)]
x0 = [0, -105AU, 0.0] # [x(0), y(0), z(0)]

for turn_param in turnparams

    rot_mat = [ 1 0 0
                0 cos(turn_param) -sin(turn_param)
                0 sin(turn_param)  cos(turn_param)]
    v0 = rot_mat*v00
    startlocvels = [vcat(x0, v0)]


for turndistoutside in turndistoutsidem
    local globalmax=globalmin=0

function lorenza(q,m,v,y_travel)

    if y_travel > -100AU #todo:drehung innen abschalten
        
        v_sol_turn=100km/s.*[0,0,1]
        
    else
        if y_travel > (-100AU-turndistoutside)
            y_travel_after100=-y_travel-100AU

            if v0[3]<0#turn upwards or dornwards
                v_sol_turn=(y_travel_after100/turndistoutside)*v0+(1-y_travel_after100/turndistoutside)*[0,0,-26km/s]
            else
                v_sol_turn=(y_travel_after100/turndistoutside)*v0+(1-y_travel_after100/turndistoutside)*[0,0,26km/s]
            end
            
        else
            v_sol_turn=v0
        end
        
    end
    
    return -q/m*cross((v-v_sol_turn),B_field(y_travel))
end

function arowdirection(y_travel)

    if y_travel > -100AU #todo:drehung innen abschalten
        
        v_sol_turn=100km/s.*[0,0,1]
        
    else
        if y_travel > (-100AU-turndistoutside)
            y_travel_after100=-y_travel-100AU

            if v0[3]<0#turn upwards or dornwards
                v_sol_turn=(y_travel_after100/turndistoutside)*v0+(1-y_travel_after100/turndistoutside)*[0,0,-26km/s]
            else
                v_sol_turn=(y_travel_after100/turndistoutside)*v0+(1-y_travel_after100/turndistoutside)*[0,0,26km/s]
            end
            
        else
            v_sol_turn=v0
        end
        
    end
    
    return v_sol_turn
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

for i in eachindex(r_dust_paramsm)
    r_dust_params=r_dust_paramsm[i]
    m_dust_params=m_dust_paramsm[i,:]



fig = Figure(resolution = (1.1*1920,1.1*1080))
#ax = Axis3(fig[1, 1], xlabel = "AU", ylabel = "AU", zlabel = "AU" ,title = "3d")
ax2 = Axis(fig[1, 1], xlabel = "AU", ylabel = "AU", title = "2d HP trajectories with incoming angle $(turn_param/pi)pi and trundistance $(turndistoutside/AU)AU")

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
    local tspan = (0, 1000day)
    

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
    lines!(ax2,label= "m=$(m)kg ,r=$(r_dust)m ,v_y=$(u0[5]/km)km/s ,v_z=$(u0[6]/km)km/s , integrtaiontime=$((tspan[2]-tspan[1])/day)days, q/m=$(round(q/m,digits=3))c/kg, turnparam_$(random_turn_param/2pi)*2pi",  Float32.(sol[2,:])/AU, Float32.(sol[3,:])/AU, linewidth = 1)
    
end
end
end
end

linex = LinRange(-100, -100, 10)
liney = LinRange(globalmin, globalmax, 10)

lines!(ax2, linex, liney, label="TS",linewidth = 2)
arrows!(ax2,[startlocvels[1][2]/AU],[startlocvels[1][3]/AU],[startlocvels[1][5]/100km],[startlocvels[1][6]/100km],normalize=true,color = :red)

Xs = LinRange(-105AU, -99AU, 15)
Ys = LinRange(0.1AU, 0.3AU, 15)                           #LinRange((globalmin+1)*AU, globalmax*AU, 10)

us = [arowdirection(x)[2] for x in Xs, y in Ys]
vs = [arowdirection(x)[3] for x in Xs, y in Ys]
strength = vec(sqrt.(us .^ 2 .+ vs .^ 2))


arrows!(Xs/AU, Ys/AU, us/AU, vs/AU, arrowsize = 10,normalize=true,color = :blue)#, lengthscale = 0.5,normalize=true, transparency= true



axislegend()
fig
save("C:\\Users\\lucac\\dustproject_clone\\DustProject\\HPplots2\\HP_startloc_$(x0/AU)_size_$(r_dust_params)_turnparam_$(turn_param/pi)pi_trundistance_$(turndistoutside/AU).png",fig)
GC.gc()
end
end
end