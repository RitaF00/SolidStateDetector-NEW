cd("/home/ritaferi/Phd/SSDdev/")

using Pkg
Pkg.activate(".")
Pkg.instantiate()

using SolidStateDetectors
using Unitful
using Plots
using LegendHDF5IO
using Printf


gr()

refinement_limits = Float32.([0.2, 0.1, 0.05, 0.02])
sim = Simulation(SSD_examples[:IVCIlayer_NODENSITY])
sim.detector = SolidStateDetector(sim.detector, contact_id=2, contact_potential=500u"V")
calculate_electric_potential!(sim,
    refinement_limits=refinement_limits,
    verbose=true,
    depletion_handling=true,
    max_tick_distance=0.1u"mm"
)

"""
pel = plot(sim.electric_potential,
    contours_equal_potential=true, full_det=true, linecolor=:white,
    levels=5,
    size=(700, 400),
    legend=false)
plot!(sim.detector, full_det=true, st=:slice, φ=0)"""

max_tick_distance = 0.1u"mm"
calculate_weighting_potential!(sim, 1,
    refinement_limits=refinement_limits,
    depletion_handling=true,
    max_tick_distance=max_tick_distance)

pwp = plot(sim.weighting_potentials[1],
    contours_equal_potential=true, linecolor=:white,
    levels=5,
    size=(700, 400),
    legend=false)
plot!(sim.detector, full_det=true, st=:slice, φ=0)
