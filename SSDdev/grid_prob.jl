cd("/home/ritaferi/Phd/SSDdev")
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using SolidStateDetectors
using Unitful
using Plots
using LegendHDF5IO



gr()


refinement_limits = [0.2, 0.1, 0.05, 0.02]
max_tick = 0.1u"mm"

sim = Simulation(SSD_examples[:InvertedCoax])
sim.detector = SolidStateDetector(sim.detector, contact_id=2, contact_potential=500u"V")

calculate_electric_potential!(sim,
    refinement_limits=refinement_limits,
    verbose=false,
    depletion_handling=true,
    grid=Grid(sim, max_tick_distance=max_tick))


calculate_weighting_potential!(sim, 1,
    refinement_limits=refinement_limits,
    verbose=false,
    depletion_handling=true,
    grid=Grid(sim, for_weighting_potential=true, max_tick_distance=max_tick))



p = plot(sim.weighting_potentials[1],
    contours_equal_potential=true,
    linecolor=:white, levels=5)
plot!(sim.detector, st=:slice, φ=0, legend=false)

savefig(p, "plots/initial_plot.png")

#=
Ep = plot(
    plot(sim.electric_potential),
    plot(sim.point_types),
    plot(sim.imp_scale),
    layout=(1, 3), size=(1200, 600))
=#
