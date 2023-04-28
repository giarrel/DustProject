include("KeptoCartfkt.jl")
include("comnames.jl")
include("constants.jl")

for comet_name in keys(comets)
    local times = range(0, stop=comets[comet_name].period, step=0.01 * day)

    # Generate orbit points
    local positions = [keplerian_to_cartesian(comet_name, times[1], time)[1] for time in times]

    # Extract X, Y, and Z components of positions
    local x_positions = [pos[1] for pos in positions]
    local y_positions = [pos[2] for pos in positions]
    local z_positions = [pos[3] for pos in positions]

    # Function to calculate distance from Sun
    distance_from_sun(x, y, z) = sqrt(x^2 + y^2 + z^2)

    # Calculate distances for all positions
    local distances = [distance_from_sun(x, y, z) for (x, y, z) in zip(x_positions, y_positions, z_positions)]

    # Find aphelion and perihelion distances
    local aphelion = maximum(distances)
    local perihelion = minimum(distances)
    local perisbdb = comets[comet_name].Q * (1 - comets[comet_name].e) / (1 + comets[comet_name].e)

    println(comet_name," (with stepsize 0.01 day):")
    println("Aphelion distance (from conversion vs sbdb): ", aphelion / AU, " AU vs ", comets[comet_name].Q / AU, " AU with rel error ", (aphelion - comets[comet_name].Q) / comets[comet_name].Q, "%")
    println("Perihelion distance distance (from conversion vs sbdb): ", perihelion / AU, " AU vs ", perisbdb / AU, " AU with rel error ", (perihelion - perisbdb) / perisbdb, "%")
    println("------------------")
end