using Plots, LinearAlgebra
include("constants.jl")
Bfield_strength=1
Bfield_strength_after=1.5e-10
Bfield_strength_before=0.5e-10

function B_field_jump(z)
    
    #B=-(v_sol_cs-v_sol)*tanh(100*z)/2-(v_sol_cs-v_sol)/2+v_sol_cs
    B=-(Bfield_strength_after-Bfield_strength_before)*tanh(100*z)/2-(Bfield_strength_after-Bfield_strength_before)/2+Bfield_strength_after

    return B
end

pts=collect(-10:0.1:10)
Plots.plot(pts,B_field_jump.(pts))

normalize([1,v_sol*sol_day/80AU,0])
v_sol*sol_day/80AU/3

