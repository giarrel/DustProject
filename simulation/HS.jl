
include("constants.jl")


d = 1/2*sol_day*v_sol_cs #distance between 2 current sheetsurfaces
Bfield_strength=1


function B_field(r,t)

    r=r+t*v_sol_cs #Field and CS moving with Solar wind

    if iseven(div(r,d)) == true #bestimmt vorzeichen von magnetfeld
        B=Bfield_strength
    else
        B=-Bfield_strength
    end

    return [0,B,0]
end

