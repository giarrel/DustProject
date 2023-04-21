using OrdinaryDiffEq
using LinearAlgebra


function launch_from_comet(initial_vel, initial_pos, tspan, beta)

    local mu = 1.327e20::Float64                            # G * M_sol [m^3 / s^2]

    # Define the second order ODE for the comet's orbit
    function orbit_ode!(ddu,du, u, p, t)
        
        beta,mu = p
        r = norm(u)
        ddu .= -mu * u / r^3 *(1-beta)

    end

    # Initial conditions and problem setup
    p = (beta,mu)

    prob = SecondOrderODEProblem(orbit_ode!, initial_vel, initial_pos, tspan, p)
    sol = solve(prob, VelocityVerlet(),dt=1000)

    x_vel_numeric = [point[1,1] for point in sol.u]
    y_vel_numeric = [point[1,2] for point in sol.u]
    z_vel_numeric = [point[1,3] for point in sol.u]

    x_pos_numeric = [point[2,1] for point in sol.u]
    y_pos_numeric = [point[2,2] for point in sol.u]
    z_pos_numeric = [point[2,3] for point in sol.u]

    return [x_vel_numeric,y_vel_numeric,z_vel_numeric] , [x_pos_numeric,y_pos_numeric,z_pos_numeric]
end