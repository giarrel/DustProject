using DifferentialEquations, CairoMakie, LinearAlgebra

#constants
const GM = 1.327e20::Float64                            # G * M_sol [m^3 / s^2]
const AU = 1.496e11::Float64                            # 1AU [m]
const yr = 3.154e7::Float64                             # 1yr [s]
const day = 24.0 * 60.0 * 60.0::Float64                 # 1d [s]

# Parameter definieren
q = 1.0
B = 5e-9
m = 1.0
w = q*B/m

# Anfangsbedingungen
v0 = [30e3, 0.0, 30e3] # [vx(0), vy(0), vz(0)]
x0 = [0.0, 0.0, 0.0] # [x(0), y(0), z(0)]
u0 = vcat(x0, v0)

# Magnetfeld
B_field = [0, 0, B]

function lorenzf(q,m,v,B_field)
    return -q/m*cross(v,B_field)
end
# Definition der Differentialgleichungen
function charged_particle!(du, u, p, t)
    du[1:3] = u[4:6]
    du[4:6] = lorenzf(q,m,u[4:6],B_field)
end

# Zeitbereich
tspan = (0.0, 15000yr)

# Analytische Lösung
x_analytic(t) = x0[1] .+ v0[1]/w*sin.(w*t) .- v0[2]/w*(1 .- cos.(w*t))
y_analytic(t) = x0[2] .+ v0[2]/w*sin.(w*t) .+ v0[1]/w*(1 .- cos.(w*t))
z_analytic(t) = x0[3] .+ v0[3]*t


# Liste der Methoden
methods_list = [RK4()]

# Anzahl der letzten Zeitschritte, die berücksichtigt werden sollen
last_N_timesteps_outside = 500000

# Schrittgröße
stepsize = 0.1day

# Wörterbuch zur Speicherung der Ergebnisse
results_dict = Dict{String, Tuple{Vector{Float64}, Vector{Float64}, Float64}}()

# Erstellen Sie einen String, um die Namen der verwendeten Methoden zu speichern
methods_used = ""

# Problem
prob = ODEProblem(charged_particle!, u0, tspan)

# Lösung für jede Methode
for method in methods_list
    local comp_time = @elapsed begin
        sol = solve(prob, method, dt=stepsize, adaptive=false)
        error_x = abs.((x_analytic(sol.t) .- sol[1,:]) ./ x_analytic(sol.t))
        error_y = abs.((y_analytic(sol.t) .- sol[2,:]) ./ y_analytic(sol.t))
        error_z = abs.((z_analytic(sol.t) .- sol[3,:]) ./ z_analytic(sol.t))
        error_magnitude = sqrt.(error_x.^2 .+ error_y.^2 .+ error_z.^2)
    end

    # Nur die letzten N Zeitschritte berücksichtigen
    if last_N_timesteps_outside != 0
        if length(sol.t) < last_N_timesteps_outside
            println("skipped ", string(method))
            continue
        end
    end

    last_N_timesteps = last_N_timesteps_outside == 0 ? length(sol.t) : last_N_timesteps_outside
    method_label = split(split(string(method), '{')[1], '(')[1]
    results_dict[method_label] = (sol.t[end-last_N_timesteps+1:end], error_magnitude[end-last_N_timesteps+1:end], comp_time)

    # Fügen Sie den Namen der aktuellen Methode der Zeichenkette hinzu
    global methods_used *= method_label * "_"
end

# Entfernen Sie das abschließende Unterstrich aus der Zeichenkette methods_used
methods_used = chop(methods_used)

# Fehler für jeden Zeitschritt plotten
f = Figure()
ax = Axis(f[1, 1], yscale=log10, title = "magnetic error, stepsize=$(stepsize/day) days", xlabel="Time [y]", ylabel="Error relative [%]")

for (method_label, (times, error, comp_time)) in results_dict
    # Filter out the values that are too small for the log plot
    valid_indices = filter(i -> error[i] > 0, 1:length(error))

    # Use only valid indices for plotting
    filtered_error = error[valid_indices]
    filtered_times = times[valid_indices] /yr

    lines!(ax, filtered_times, filtered_error, label="$(method_label) ($(comp_time) s)")
end

axislegend(ax)

# Save the plot
save("error_plots/magnetic_error_stepsize$(stepsize/day)_days_inttime_$(tspan[2]/yr)yr_methods_$(methods_used).png",f)