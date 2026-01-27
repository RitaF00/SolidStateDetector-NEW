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

#-------- ICPC vuoto -----------
#save_sim_path = "saved_simulation/sim.h5"
#-------- ICPC ILM -----------
save_sim_path = "saved_simulation/sim_ILM.h5"



# ------------------------
# SCELTA 
# ------------------------
# true  → ricalcola sempre e sovrascrivi
# false → usa il file salvato (se esiste)
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

#------------ qui inizia la simulazione del weighting potential -------
# inizializzo lo stato con una griglia più grande da 1 mm
new_grid = Grid(sim, max_tick_distance=1u"mm")
apply_initial_state!(sim, WeightingPotential, 1, new_grid) # optional
SolidStateDetectors.update_till_convergence!(sim, WeightingPotential, 1, depletion_handling=true, verbose=true)

# inizio con il refinement della griglia costruita con 1 mm
refinement_limits = Float32.([0.2, 0.1, 0.5, 0.2])
n_refinement_steps = length(refinement_limits)
for iref in 1:n_refinement_steps
    println("Refinement step $(refinement_limits[iref])")
    max_diff_array = (refinement_limits[iref], refinement_limits[iref], refinement_limits[iref])
    SolidStateDetectors.refine!(sim, WeightingPotential, 1, max_diff_array)
    SolidStateDetectors.update_till_convergence!(sim, WeightingPotential, 1, depletion_handling=true, verbose=true)
end

d = 0.5u"mm"
sim.weighting_potentials[1] =
    SolidStateDetectors.refine_scalar_potential_by_tick_distance(
        sim.weighting_potentials[1],
        (d, d, d),
        (0.0001u"mm", 0.0001u"mm", 0.0001u"mm")
    );