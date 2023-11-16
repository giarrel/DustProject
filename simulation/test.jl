using DifferentialEquations, GLMakie, LinearAlgebra
include("constants.jl")
include("duststartHS.jl")


# Anfangsbedingungen
v0 = [0,26,0]* km/s # [vx(0), vy(0), vz(0)] TODO: RICHTUNG WIE DAS TEILCHEN REINFLIEGT
x0 = [0, -101AU, 0.0] # [x(0), y(0), z(0)]
startlocvels = [vcat(x0, v0)]

turndistoutside=1AU

xs = [-100AU]
ys = [0]
linex = LinRange(-100AU, -100AU, 16)
liney = LinRange(-1.1AU, 1.1AU, 16)

function arowdirection(y_travel)
    phi=pi/2
    if y_travel > -100AU #todo:drehung innen abschalten
        
        v_sol_turn=100km/s.*[0,0,1]
        
    else
        if y_travel > (-100AU-turndistoutside)
            y_travel_after100=-y_travel-100AU
            
            v_sol_turn=(y_travel_after100/turndistoutside)*v0+(1-y_travel_after100/turndistoutside)*[0,0,26km/s]
        else
            v_sol_turn=v0
        end
        
    end
    
    return v_sol_turn
end

us = [arowdirection(x)[2] for x in xs, y in ys]
vs = [arowdirection(x)[3] for x in xs, y in ys]


fig = Figure(resolution = (1.1*1920,1.1*1080))

ax2 = Axis(fig[1, 2], xlabel = "AU", ylabel = "AU", title = "Gasflow 2d with turn over distance of $(turndistoutside/AU)AU")

arrows!(xs/AU, ys/AU, us/AU, vs/AU, arrowsize = 10, lengthscale = 0.1,
     normalize=true)

lines!(linex/AU,liney/AU,label="Shock")
axislegend()
fig
#save("C:\\Users\\lucac\\dustproject_clone\\DustProject\\HPplots2\\HP_startloc_$(x0/AU)_size_$(r_dust_params)_turnparam_$(random_turn_param/2pi)_trumpara_$(turndistoutside/AU).png",fig)
GC.gc()

sin(500pi)