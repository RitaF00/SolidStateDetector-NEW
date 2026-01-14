cd("/home/ritaferi/Phd/SSDdev/")
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using SolidStateDetectors
using Unitful
using Plots
using LegendHDF5IO

gr()

refinement_limits = [0.2, 0.1, 0.05, 0.02]

#-------- ICPC vuoto -----------
save_sim_path = "saved_simulation/sim.h5"
sim = ssd_read(save_sim_path, Simulation)

max_tick_distance = 0.5u"mm"

calculate_weighting_potential!(sim, 1,
    refinement_limits=refinement_limits,
    depletion_handling=true,
    grid=Grid(sim,
        for_weighting_potential=true,
        max_tick_distance=max_tick_distance))



x_min = 0.015
x_max = 0.020
y_min = 0.02
y_max = 0.04
x_rect = [x_min, x_max, x_max, x_min, x_min]
y_rect = [y_min, y_min, y_max, y_max, y_min]

plot!(plt[i], sim.weighting_potentials[1],
    contours_equal_potential=true,
    linecolor=:white,
    levels=5,
    title="max tick distance = $max_tick_distance")

# Plot del detector
plot!(plt[i], sim.detector, st=:slice, φ=0, legend=false)

# Rettangolo
plot!(plt[i], x_rect, y_rect,
    seriestype=:shape,
    linecolor=:white,
    lw=1.5,
    fillalpha=0,
    label="")

display(plt)
savefig(plt, "solo_rett.png")