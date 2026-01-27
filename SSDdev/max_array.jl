cd("/home/ritaferi/Phd/SSDdev/")
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using SolidStateDetectors
using Unitful
using Plots

gr()

sim = Simulation(SSD_examples[:IVCIlayer])
sim.detector = SolidStateDetector(sim.detector, contact_id=2, contact_potential=3500u"V")

calculate_electric_potential!(sim,
    refinement_limits=[0.2],
    verbose=true, #  boolean in the output is produced or not
    depletion_handling=true)