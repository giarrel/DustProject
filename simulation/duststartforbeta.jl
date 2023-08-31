using LinearAlgebra

include("constants.jl")


function startparams_beta(no_of_starts,vd,beta)

    
    l,b  = (deg2rad(0),deg2rad(0))
    x,y,z = cos(b)*cos(l),cos(b)*sin(l),sin(b)
    start_locvel=zeros(6,no_of_starts)
    for i in 1:no_of_starts
        start_locvel[1:3,i]= [0AU,0AU,0AU]+80AU*[x,y,z]+2AU*(rand()-0.5)*[y,-x,0]/norm([y,-x,0])
        v=sqrt(vd^2+2GM*(1-beta)/norm(start_locvel[1:3,i]))
        start_locvel[4:6,i]=v*[-x,-y,-z]
    end
    return start_locvel
end