using SolidStateDetectors, Unitful

sim = Simulation{Float64}(SSD_examples[:BEGe])
sim.detector = SolidStateDetector(sim.detector, contact_id=2, contact_potential=3500u"V")
calculate_electric_potential!(sim)
calculate_electric_field!(sim)

for i in 1:2
    calculate_weighting_potential!(sim, i)
end


locations = [CartesianPoint(0.035, 0, 0.02)]
energies = [1000u"keV"]
evt = Event(locations, energies)
simulate!(evt, sim)
