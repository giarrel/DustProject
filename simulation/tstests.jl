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

B_field(-81AU)
B_field(-79AU)


function lorenza(y_travel)
    y_travel_v_jump=y_travel+80AU
    return v_sol_jump = -(v_sol_cs-v_sol)*tanh(100*y_travel_v_jump)/2-(v_sol_cs-v_sol)/2+v_sol_cs
    #return -q/m*cross((v-[0,v_sol_jump,0]),B_field(y_travel))
end

lorenza(-81AU)
lorenza(-79AU)