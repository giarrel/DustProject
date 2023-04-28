include("KeptoCartfkt.jl")
include("comnames.jl")
include("constants.jl")

for comet_name in keys(comets)
    local period = comets[comet_name].period
    local times_all = range(0, stop=100 * period, step=0.1 * day)
    local timebegin_last = (100 - 1) * period  # Start at the second to last period
    local timeend_last = 20 * period          # End at the last period

    # Generate orbit points for all times
    local positions_all = [keplerian_to_cartesian(comet_name, times_all[1], time)[1] for time in times_all]

    # Find the indices corresponding to the last period
    idx_last_period = findall(t -> t >= timebegin_last, times_all)

    # Extract positions for the last period
    local positions_last = positions_all[idx_last_period]

    # Extract X, Y, and Z components of positions
    local x_positions = [pos[1] for pos in positions_last]
    local y_positions = [pos[2] for pos in positions_last]
    local z_positions = [pos[3] for pos in positions_last]

    # Function to calculate distance from Sun
    distance_from_sun(x, y, z) = sqrt(x^2 + y^2 + z^2)

    # Calculate distances for all positions
    local distances = [distance_from_sun(x, y, z) for (x, y, z) in zip(x_positions, y_positions, z_positions)]

    # Find aphelion and perihelion distances
    local aphelion = maximum(distances)
    local perihelion = minimum(distances)
    local perisbdb = comets[comet_name].Q * (1 - comets[comet_name].e) / (1 + comets[comet_name].e)

    println(comet_name," (with stepsize 0.1 day):")
    println("Aphelion distance (from conversion vs sbdb): ", aphelion / AU, " AU vs ", comets[comet_name].Q / AU, " AU with rel error ", (aphelion - comets[comet_name].Q) / comets[comet_name].Q, "%")
    println("Perihelion distance distance (from conversion vs sbdb): ", perihelion / AU, " AU vs ", perisbdb / AU, " AU with rel error ", (perihelion - perisbdb) / perisbdb, "%")
    println("------------------")
end
println("-----------------")