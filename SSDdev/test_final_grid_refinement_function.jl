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


save_sim_path = "saved_simulation/sim_ILM.h5"

recalculate = false   # <<<<<<<< CAMBIA QUI
if isfile(save_sim_path) && !recalculate
    println("⚡ Upoload simulation saved in : $save_sim_path")
    sim = ssd_read(save_sim_path, Simulation)

else
    println("🔧 New simulation for the electric potential...")
    max_tick_distance = 0.5u"mm"
    #sim = Simulation(SSD_examples[:InvertedCoax])
    sim = Simulation(SSD_examples[:IVCIlayer])
    sim.detector = SolidStateDetector(sim.detector, contact_id=2, contact_potential=500u"V")
    grid = Grid(sim, max_tick_distance=max_tick_distance)


    calculate_electric_potential!(sim,
        refinement_limits=[0.2, 0.1, 0.05, 0.02],
        verbose=true, #  boolean in the output is produced or not
        depletion_handling=true,  # motiplica epsilon_r per un fattore f nelle regioni non svuotate
        grid=grid)

    println("💾 Saving simulation in $save_sim_path")
    ssd_write(save_sim_path, sim)
end



#max_tick_distance = (0.5u"mm", 1u"rad", 0.5u"mm")
max_tick_distance = 0.5u"mm"

refinement_limits = [0.2, 0.1, 0.05, 0.02]
# inizio del time

SolidStateDetectors._calculate_potential_max_tick_refinement!(sim,
    WeightingPotential,
    1,
    max_tick_distance=max_tick_distance,
    depletion_handling=true,
    verbose=true,
    refinement_limits=refinement_limits)


# fine del tempo
t1 = time_ns()
elapsed = (t1 - t0) / 1e9  # secondi


println("========== Total run time = $elapsed s ============")


println("calcolo del minimo del potenziale")
wp = sim.weighting_potentials[1]
r_ticks = wp.grid.axes[1].ticks
z_ticks = wp.grid.axes[3].ticks

rmin, rmax = ustrip.(u"m", (0.015u"m", 0.02u"m"))
zmin, zmax = ustrip.(u"m", (0.02u"m", 0.03u"m"))

r_inds = findall(r -> rmin ≤ r ≤ rmax, r_ticks)
z_inds = findall(z -> zmin ≤ z ≤ zmax, z_ticks)

isempty(r_inds) && error("No points found in the r interval")
isempty(z_inds) && error("No points found in the z interval")

vals = wp.data[r_inds, 1, z_inds]
min_wp = minimum(vals)
println(">>> Min WeightingPotential in the test volume = $min_wp")