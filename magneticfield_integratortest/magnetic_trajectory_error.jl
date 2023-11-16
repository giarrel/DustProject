using DifferentialEquations, CairoMakie, LinearAlgebra

#constants
const GM = 1.327e20::Float64                            # G * M_sol [m^3 / s^2]
const AU = 1.496e11::Float64                            # 1AU [m]
const yr = 3.154e7::Float64                             # 1yr [s]
const day = 24.0 * 60.0 * 60.0::Float64                 # 1d [s]

# Parameter definieren
q = 1.0
B = 1e-5
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
tspan = (0.0, 200yr)

# Analytische Lösung https://suli.pppl.gov/2018/course/Fox_SULI_2018.pdf

x_analytic(t) = v0[1]/w*sin.(w*t)
y_analytic(t) = v0[1]/w*(1 .- cos.(w*t))
z_analytic(t) = v0[3]*t

#=
x_analytic(t) = x0[1] .+ v0[1]/w*sin.(w*t) .- v0[2]/w*(1 .- cos.(w*t))
y_analytic(t) = x0[2] .+ v0[2]/w*sin.(w*t) .+ v0[1]/w*(1 .- cos.(w*t))
z_analytic(t) = x0[3] .+ v0[3]*t
=#

# Liste der Metho
methods_list = [
#=
Euler(),

Midpoint(),
ImplicitEuler(),
ImplicitMidpoint(),
Trapezoid(),
RK4(),
Vern9(),
AutoVern9(Rodas5P())=#

    #=Adams/multistep methods
    AB3(),
    AB4(),
    AB5(),
    ABM32(),
    ABM43(),=#
    ABM54()
#

    ]
    

# Anzahl der letzten Zeitschritte, die berücksichtigt werden sollen
last_N_timesteps_outside = 00000

# Schrittgröße
stepsize = 0.05day
adap=false
reltol=1e-13

# Wörterbuch zur Speicherung der Ergebnisse
results_dict = Dict{String, Tuple{Vector{Float64}, Vector{Float64}, Float64}}()

# Erstellen Sie einen String, um die Namen der verwendeten Methoden zu speichern
methods_used = ""

# Problem
prob = ODEProblem(charged_particle!, u0, tspan)

# Lösung für jede Methode
for method in methods_list
    local comp_time = @elapsed begin
        sol = solve(prob, method, dt=stepsize, adaptive=adap,reltol=reltol)
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
if adap == false
    ax = Axis(f[1, 1], yscale=log10, title = "magnetic error, step size = $(stepsize/day) days", xlabel="Time [y]", ylabel="Error relative [%]")
else
    ax = Axis(f[1, 1], yscale=log10, title = "magnetic error, reltol = $(reltol)", xlabel="Time [y]", ylabel="Error relative [%]")
end


for (method_label, (times, error, comp_time)) in results_dict
    # Filter out the values that are too small for the log plot
    valid_indices = filter(i -> error[i] > 1e-18, 1:length(error))

    # Use only valid indices for plotting
    filtered_error = error[valid_indices]
    filtered_times = times[valid_indices] /yr

    lines!(ax, filtered_times, filtered_error, label="$(method_label) ($(comp_time) s)")
end

axislegend(ax)

# Save the plot
if adap == false
    save("DustProject/error_plots_new/magnetic_error_B$(B)_stepsize$(stepsize/day)_days_inttime_$(tspan[2]/yr)yr_methods_$(methods_used).pdf",f)
else
    save("DustProject/error_plots_ada2/magnetic_error_B$(B)_reltol$(reltol)_days_inttime_$(tspan[2]/yr)yr_methods_$(methods_used).pdf",f)
end