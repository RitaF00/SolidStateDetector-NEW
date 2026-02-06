
cd("/home/ritaferi/Phd/SSDdev/")
using Pkg
Pkg.activate(".")
Pkg.instantiate()



using SolidStateDetectors
using Geant4
using Plots
using Unitful
using StatsBase
using LegendHDF5IO
T = Float32

sim = Simulation(SSD_examples[:IVCIlayer])
calculate_electric_potential!(sim)
calculate_electric_field!(sim)
for i in 1:2
    calculate_weighting_potential!(sim, i)
end

println("🔆 Simulating an event...")
evt = Event([CartesianPoint{T}(0.02, 0, 0.05)], [1000u"keV"])
simulate!(evt, sim)

p = plot(u"mm", u"mm", u"mm")
plot!(sim.detector, st=:slice, φ=0u"°", size=(500, 500), label="")
plot!(evt.drift_paths, linewidth=2, linestyle=:dash, markersize=5)

savefig(p, "IVC_event.png")
