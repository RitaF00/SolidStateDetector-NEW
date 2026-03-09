cd("/home/ritaferi/Phd/SSDdev/")
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using SolidStateDetectors
using Unitful
using Plots
using LegendHDF5IO
using HDF5

gr()



sim = Simulation(SSD_examples[:InfiniteCoaxialCapacitor])
sim.detector = SolidStateDetector(sim.detector, contact_id=2, contact_potential=3500u"V")
plot(
    plot(sim.detector, seriestype=:samplesurface, title=":samplesurface",
        markersize=1, n_samples=50),
    size=(400, 400), legend=false, ticks=false,
    guide="", zlims=(-0.005, 0.1), axis=false
)