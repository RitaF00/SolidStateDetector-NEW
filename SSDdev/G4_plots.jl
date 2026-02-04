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

println("ESCI E RIENTRA ONGI VOLTA")

sim = Simulation{T}(SSD_examples[:InvertedCoaxInCryostat])

source_1 = MonoenergeticSource(
    "gamma",                              # Type of particle beam
    2.615u"MeV",                          # Energy of particle
    CartesianPoint(0.065, 0., 0.05),      # Location of the source
    CartesianVector(-1, 0, 0),              # Direction of the source
    10u"°"                                # Opening angle of the source emission
)

app = G4JLApplication(sim, source_1, verbose=false);
N_events = 50000
events = run_geant4_simulation(app, N_events)

p = plot(sim.detector, show_passives=false, size=(500, 500), fmt=:png)
plot!(source_1)
pts = [
    CartesianPoint(p.x, p.y, p.z)
    for p in events[1:1000].pos.data
]

plot!(pts, ms=0.5, msw=0, color=:black, label="")
savefig(p, "InvertedCoax_G4_events.png")


h = fit(Histogram, ustrip.(u"keV", sum.(events.edep)), Weights(fill(10, length(events.edep))), 0:10:3000)


p = plot(
    h,
    bins=500,
    yscale=:log10,
    st=:step,
    label="Simulated data",
    framestyle=:box,
    grid=true,
    gridalpha=0.3,
    tickfont=font(14),                    # dimensione numeri asse
    xlabel="E (keV)",
    ylabel="counts",
    guidefontsize=18,
    legendfontsize=14,                    # dimensione testo legenda
    titlefontsize=20,                     # dimensione del titolo (se presente)
    lw=3                                  # spessore istogramma
)

xlims!(0, 3000)
xticks!(0:500:3000)

savefig(p, "simulated_spectrum.png")