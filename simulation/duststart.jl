using LinearAlgebra

include("constants.jl")


function startparams(no_of_starts,v=26,l=255.41,b=5.03)

    l,b  = (deg2rad(l),deg2rad(b))
    x,y,z = cos(b)*cos(l),cos(b)*sin(l),sin(b)

    start_locvel=zeros(6,no_of_starts)
    for i in 1:no_of_starts
        start_locvel[4:6,i]=v*[-x,-y,-z]  #start_locvel[4:6,i]=(v+(0.5-rand()))*[-x,-y,-z]*1e3
        start_locvel[1:3,i]= [0AU,0AU,0AU]+80AU*[x,y,z]+1AU*(rand()-0.5)*[y,-x,0]/norm([y,-x,0])#+80AU*(rand()-0.5)*cross([x,y,z],[y,-x,0]/norm([y,-x,0]))
    end
    return start_locvel
end