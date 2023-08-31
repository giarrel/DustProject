using DifferentialEquations, GLMakie, LinearAlgebra,NLsolve, Plots
include("constants.jl")

d = 1/2*sol_day*v_sol_cs #distance between 2 current sheetsurfaces
Bfield_strength=1.5e-10
v_dust_init=26*km/s
r_dust=1e-9
m_dust=1e-22
q=5*4*pi*eps_0*r_dust

solve_for_zero(t) = d+m_dust*(v_dust_init-v_sol_cs)/(abs(q)*Bfield_strength)*sin(abs(q)*Bfield_strength*t[1]/m_dust)#-m_dust*v_sol_cs/(abs(q)*Bfield_strength)*sin(abs(q)*Bfield_strength*t[1]/m_dust)#v_sol_cs*t[1]
solve_for_zero2(t) = m_dust*v_dust_init/(abs(q)*Bfield_strength)*sin(abs(q)*Bfield_strength*t[1]/m_dust)*(abs(q)*Bfield_strength/m_dust)^2


nlsolve(solve_for_zero,[955065.6000001743])


t=collect(1:1000:10000000)

Plots.plot(t,solve_for_zero.(t))
