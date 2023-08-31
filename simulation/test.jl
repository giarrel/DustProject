using Plots
Bfield_strength=1

function B_field_min(z)
    

    B=-2*Bfield_strength*tanh(100*z)/2

    return B
end

pts=collect(-10:0.1:10)
Plots.plot(pts,B_field_min.(pts))