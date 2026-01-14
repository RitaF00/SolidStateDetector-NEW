

using Pkg
Pkg.activate(".")
Pkg.instantiate()

using SolidStateDetectors
using Unitful
using Plots
using LegendHDF5IO

gr()

sim = Simulation(SSD_examples[:InvertedCoax])
sim.detector = SolidStateDetector(sim.detector, contact_id=2, contact_potential=500u"V")

calculate_electric_potential!(sim, refinement_limits=[0.2, 0.1, 0.05, 0.01], verbose=false, depletion_handling=true, grid=Grid(sim, max_tick_distance=0.5u"mm"))

calculate_weighting_potential!(sim, 1, refinement_limits=[0.2, 0.1, 0.05, 0.01], verbose=false, depletion_handling=true, grid=Grid(sim, max_tick_distance=0.5u"mm"))

p = plot(sim.weighting_potentials[1],
    contours_equal_potential=true,
    linecolor=:white, levels=5)
plot!(sim.detector, st=:slice, φ=0, legend=false)