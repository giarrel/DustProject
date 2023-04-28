#https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf
include("constants.jl")
#convert a traditional set of Keplerian Orbit Elements
#– Semi-major axis a [m]
#– Eccentricity e [1]
#– Argument of periapsis ω [rad]
#– Longitude of ascending node (LAN) Ω [rad]
#– Inclination i [rad]
#– Mean anomaly M0 = M(t0) [rad] at epoch t0 [s]
#at epoch t
#to position vector r(t) [m] and r_dot(t) [m/s]


function keplerian_to_cartesian(comet_name::String,t0,t; tol=1e-10)
    
    comet = comets[comet_name]
    a, e, ω, Ω, i, M0 = comet.a, comet.e, comet.peri, comet.node, comet.i, comet.M0
    
    t=(t-t0)

    function M_t()
        return (M0 + t * sqrt(GM / a^3))%(2pi)
    end

    function E_t(tolerance)
        E_i = M_t()
        while abs(E_i - e * sin(E_i) - M_t()) > tolerance
            E_i -= (E_i - e * sin(E_i) - M_t()) / (1 - e * cos(E_i))
        end
        return E_i
    end

    E = E_t(tol)

    function true_anomaly()
        arg1 = sqrt(1 + e) * sin(E / 2)
        arg2 = sqrt(1 - e) * cos(E / 2)
        return 2 * atan(arg1, arg2)
    end

    ν = true_anomaly()

    function r_c()
        return a * (1 - e * cos(E))
    end

    rc = r_c()

    function orbital_vectors()
        # Position vector o(t)
        ox = rc * cos(ν)
        oy = rc * sin(ν)
        oz = 0.0
        o = [ox, oy, oz]
    
        # Velocity vector o_dot(t)
        o_dot_x = sqrt(GM * a) / rc * (-sin(E))
        o_dot_y = sqrt(GM * a) / rc * (sqrt(1 - e^2) * cos(E))
        o_dot_z = 0.0
        o_dot = [o_dot_x, o_dot_y, o_dot_z]
    
        return o, o_dot
    end

    o, o_dot = orbital_vectors()

    function orbital_to_inertial()
        ox, oy = o[1], o[2]
        o_dot_x, o_dot_y = o_dot[1], o_dot[2]
    
        # Calculate r(t)
        r_x = ox * (cos(ω) * cos(Ω) - sin(ω) * cos(i) * sin(Ω)) - oy * (sin(ω) * cos(Ω) + cos(ω) * cos(i) * sin(Ω))
        r_y = ox * (cos(ω) * sin(Ω) + sin(ω) * cos(i) * cos(Ω)) + oy * (cos(ω) * cos(i) * cos(Ω) - sin(ω) * sin(Ω))
        r_z = ox * (sin(ω) * sin(i)) + oy * (cos(ω) * sin(i))
    
        # Calculate r_dot(t)
        r_dot_x = o_dot_x * (cos(ω) * cos(Ω) - sin(ω) * cos(i) * sin(Ω)) - o_dot_y * (sin(ω) * cos(Ω) + cos(ω) * cos(i) * sin(Ω))
        r_dot_y = o_dot_x * (cos(ω) * sin(Ω) + sin(ω) * cos(i) * cos(Ω)) + o_dot_y * (cos(ω) * cos(i) * cos(Ω) - sin(ω) * sin(Ω))
        r_dot_z = o_dot_x * (sin(ω) * sin(i)) + o_dot_y * (cos(ω) * sin(i))
    
        return [r_x, r_y, r_z], [r_dot_x, r_dot_y, r_dot_z]
    end
    
    r, r_dot = orbital_to_inertial()
    
    return r, r_dot
end