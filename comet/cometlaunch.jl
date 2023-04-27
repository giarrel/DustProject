using OrdinaryDiffEq
using DifferentialEquations
using LinearAlgebra


function launch_from_comet(initial_vel, initial_pos, tspan, beta, method; adaptive=false, stepsize=1e5, abstol=1e-3, reltol=1e-3)

    local mu = 1.327e20::Float64                            # G * M_sol [m^3 / s^2]

    # Define the second order ODE for the comet's orbit
    function orbit_ode(ddu,du, u, p, t)
        
        beta,mu = p
        r = norm(u)
        ddu .= -mu * u / r^3 *(1-beta)

    end

    # Initial conditions and problem setup
    p = (beta,mu)

    prob = SecondOrderODEProblem(orbit_ode, initial_vel, initial_pos, tspan, p)
    sol = solve(prob, method, adaptive=adaptive, dt=stepsize, abstol=abstol, reltol=reltol)

    x_vel_numeric = [point[1,1] for point in sol.u]
    y_vel_numeric = [point[1,2] for point in sol.u]
    z_vel_numeric = [point[1,3] for point in sol.u]

    x_pos_numeric = [point[2,1] for point in sol.u]
    y_pos_numeric = [point[2,2] for point in sol.u]
    z_pos_numeric = [point[2,3] for point in sol.u]

    times = sol.t
    
    return ([x_vel_numeric,y_vel_numeric,z_vel_numeric] , [x_pos_numeric,y_pos_numeric,z_pos_numeric], times)
end



function launch_from_comet_1ord(initial_vel, initial_pos, tspan, beta, method; adaptive=false, stepsize=1e5, abstol=1e-3, reltol=1e-3)

    local mu = 1.327e20::Float64                            # G * M_sol [m^3 / s^2]

    # Define the first order ODE system for the comet's orbit
    function orbit_ode!(du, u, p, t)
        
        beta, mu = p
        r = norm(u[1:3])
        du[1:3] = u[4:6]
        du[4:6] = -mu * u[1:3] / r^3 * (1 - beta)

    end

    # Initial conditions and problem setup
    p = (beta, mu)
    u0 = vcat(initial_pos, initial_vel)

    prob = ODEProblem(orbit_ode!, u0, tspan, p)
    sol = solve(prob, method, adaptive=adaptive, dt=stepsize, abstol=abstol, reltol=reltol)

    # Extract the results
    x_vel_numeric = [point[4] for point in sol.u]
    y_vel_numeric = [point[5] for point in sol.u]
    z_vel_numeric = [point[6] for point in sol.u]

    x_pos_numeric = [point[1] for point in sol.u]
    y_pos_numeric = [point[2] for point in sol.u]
    z_pos_numeric = [point[3] for point in sol.u]

    times = sol.t

    return ([x_vel_numeric, y_vel_numeric, z_vel_numeric], [x_pos_numeric, y_pos_numeric, z_pos_numeric], times)
end
