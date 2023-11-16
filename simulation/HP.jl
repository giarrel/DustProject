using DifferentialEquations, GLMakie, LinearAlgebra
include("constants.jl")
include("duststartHS.jl")


Bfield_strength_after=6.0e-10
Bfield_strength_before=1.5e-10

random_turn_param=0.1*2pi
rot_mat = [ cos(random_turn_param) 0 -sin(random_turn_param)
            0                      1        0
            sin(random_turn_param) 0 cos(random_turn_param)]

B_direction_inside=[1,0,0]
B_direction_outside=*(rot_mat,B_direction_inside)


function B_field(y_travel)
    y_travel_adjust=y_travel+100AU
    B=-(Bfield_strength_after-Bfield_strength_before)*tanh(100*y_travel_adjust)/2-(Bfield_strength_after-Bfield_strength_before)/2+Bfield_strength_after

    if y_travel>-100AU
        return B.*B_direction_inside
    else
        return B.*B_direction_outside
    end
end

B_field(-101AU)
B_field(-99AU)
B_field(50AU)

function lorenza(q,m,v,y_travel)
    phi=pi/2
    if y_travel > -100AU
        if y_travel < -90AU
            y_travel_after100=y_travel+100AU
            phi = y_travel_after100/10AU*pi/2
            v_sol_turn=-100km/s.*[cos(phi),sin(phi),0]
        else
            v_sol_turn=-100km/s.*[0,1,0]
        end
    else
        v_sol_turn=[0,0,0]
    end
    
    return -q/m*cross((v-v_sol_turn),B_field(y_travel))#-[0,v_sol_jump,0] relativgeschwindigkeit nur innen HP, TODO:man kann als was wäre sich anschauen was wenn ausserhalb eingefrohren beide seiten parallel zum hp aber nicht unbedingt parallel zueinander
end


beta_params = [0]
r_dust_params=[1e-6]
m_dust_params=[1e-14,1e-15,1e-16]
#r_dust_params=[1e-7]
#m_dust_params=[1e-17,1e-18,1e-19]
#r_dust_params=[1e-8]
#m_dust_params=[1e-20,1e-21,1e-22]
#r_dust_params=[1e-9]
#m_dust_params=[1e-22,1e-23,1e-24]



# Anfangsbedingungen
v0 = [0,26,0]* km/s # [vx(0), vy(0), vz(0)] TODO: RICHTUNG WIE DAS TEILCHEN REINFLIEGT
x0 = [0, -101AU, 0.0] # [x(0), y(0), z(0)]
startlocvels = [vcat(x0, v0)]



fig = Figure(resolution = (1.1*1920,1.1*1080))
ax = Axis3(fig[1, 1], xlabel = "AU", ylabel = "AU", zlabel = "AU" ,title = "HP dust trajectory with magnetic field turn ed by $(random_turn_param/pi)pi 3d")
ax2 = Axis(fig[1, 2], xlabel = "AU", ylabel = "AU", title = "HP dust trajectory with magnetic field turn ed by $(random_turn_param/pi)pi 2d")

for u0 in startlocvels
for m in m_dust_params
for r_dust in r_dust_params
for beta in beta_params

    local q=5*4*pi*eps_0*r_dust

    

    # Definition der Differentialgleichungen
    function charged_particle!(du, u, p, t)
        r = norm(u[1:3])
        du[1:3] = u[4:6]
        du[4:6] = lorenza(q,m,u[4:6],u[2])#-GM * u[1:3] / r^3 * (1 - beta)
    end

    # Zeitbereich
    local tspan = (0, 200day)
    

    # Problem und Lösung
    local prob = ODEProblem(charged_particle!, u0, tspan)
    local sol = solve(prob,Vern9(),dt=0.001*day,adaptive=true,reltol=1e-20)

    # Positionen plotten
    lines!(ax,label= string.(m), Float32.(sol[1,:])/AU, Float32.(sol[2,:])/AU, Float32.(sol[3,:])/AU, linewidth = 1)
    lines!(ax2,label= "m=$(m)kg ,r=$(r_dust)m ,v_y=$(u0[5]/km)km/s ,v_z=$(u0[6]/km)km/s , integrtaiontime=$((tspan[2]-tspan[1])/day)days, q/m=$(round(q/m,digits=3))c/kg, turnparam_$(random_turn_param/2pi)2pi",  Float32.(sol[2,:])/AU, Float32.(sol[3,:])/AU, linewidth = 1)
    
end
end
end
end

axislegend()
fig
save("C:\\Users\\lucac\\dustproject_clone\\DustProject\\HPplots\\HP_startloc_$(x0/AU)_size_$(r_dust_params)_turnparam_$(random_turn_param/pi)pi.png",fig)
GC.gc()