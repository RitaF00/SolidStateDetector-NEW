cd("/home/ritaferi/Phd/SSDdev")
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using SolidStateDetectors
using Unitful
using Plots
using LegendHDF5IO



#----------- 0.05 tick grid ------------
max_tick_distance = 0.05u"mm"
#refinement_limits = [0.2, 0.1, 0.05, 0.02]
refinement_limits = [1]
V_bias = 3500u"V"


# true  → ricalcola sempre e sovrascrivi
# false → usa il file salvato (se esiste)
recalculate = false   # <<<<<<<< CAMBIA QUI
save_sim_path = "saved_simulation/sim.05"

if isfile(save_sim_path) && !recalculate
    println("⚡ Upoload simulation saved in : $save_sim_path")
    sim05 = ssd_read(save_sim_path, Simulation)

else
    println("🔧 New simulation for the electric potential...")

    sim05 = Simulation(SSD_examples[:InvertedCoax])
    sim05.detector = SolidStateDetector(sim05.detector, contact_id=2, contact_potential=V_bias)
    grid = Grid(sim05, max_tick_distance=max_tick_distance)


    # ATTENZIONE  A VEDERE SE IL POTENZIALE E' DEPLETED O NO
    calculate_electric_potential!(sim05,
        refinement_limits=missing,
        verbose=false, #  boolean in the output is produced or not
        depletion_handling=false,  # se true : motiplica epsilon_r per un fattore f nelle regioni non svuotate
        grid=Grid(sim05, max_tick_distance=max_tick_distance))

    println("💾 Saving simulation in $save_sim_path")
    sim05 = ssd_write(save_sim_path, sim)
end


final_plot = plot(
    plot(sim05.electric_potential, φ=0), # initial electric potential (boundary conditions)
    plot(sim05.point_types), # map of different point types: fixed point / inside or outside detector volume / depleted/undepleted
    plot(sim05.imp_scale), # impurity scale
    plot(sim05.q_eff_imp), # charge density distribution
    layout=(1, 4), size=(3000, 1000),
    title="grid tick $max_tick_distance and V_bias $V_bias "
)
savefig(final_plot, "plots/electric_pot_prob/0.05mm_grid.png")
