cd("/home/ritaferi/Phd/SSDdev")
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using SolidStateDetectors
using Unitful
using Plots
using LegendHDF5IO



gr()

println("Initial plot")
refinement_limits = [0.2, 0.1, 0.05, 0.02]
max_tick = 0.1u"mm"

using SolidStateDetectors
using Unitful
using Plots

gr()

println("Initial plot")

# Inizializzo la simulazione
sim = Simulation(SSD_examples[:IVCIlayer])
sim.detector = SolidStateDetector(sim.detector, contact_id=2, contact_potential=500u"V")

# Definisco i max_tick da testare
max_ticks = 0.1u"mm"
convergence_limit = 1e-7
n_iterations_between_checks = 10000



grid = Grid(sim, max_tick_distance=max_ticks)
println("calculatong electric potential...")
# Calcolo potenziale elettrico
calculate_electric_potential!(sim;
    refinement_limits=refinement_limits,
    verbose=false,
    depletion_handling=true,
    grid=grid)

println("Calculating weighting potential")
# Calcolo potenziale pesato usando la simulazione dell'electric potential appena calcolata
calculate_weighting_potential!(sim, 1;
    n_iterations_between_checks=n_iterations_between_checks,
    convergence_limit=convergence_limit,
    refinement_limits=refinement_limits,
    verbose=false,
    depletion_handling=true,
    grid=Grid(sim, for_weighting_potential=true, max_tick_distance=max_ticks))


p = plot(sim.weighting_potentials[1],
    contours_equal_potential=true,
    linecolor=:white, levels=5,
    title="number iteration bewteen checks = $(n_iterations_between_checks) and convergence limit = $(convergence_limit)")
plot!(sim.detector, st=:slice, φ=0, legend=false)


savefig(p, "plots/ILM_10-7_10000.png")
