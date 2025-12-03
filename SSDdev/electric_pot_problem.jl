cd("/home/ritaferi/Phd/SSDdev")
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using SolidStateDetectors
using Unitful
using Plots
using LegendHDF5IO





refinement_limits = 0.2
V_bias = 3500u"V"

println("Refinement limiits $refinement_limits")

#----------- 0.05 tick grid ------------
max_tick_distance = 0.05u"mm"

println("Eletric potential for grid $max_tick_distance")
# true  → ricalcola sempre e sovrascrivi
# false → usa il file salvato (se esiste)
recalculate = false   # <<<<<<<< CAMBIA QUI
save_sim_path = "saved_simulation/sim05"

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
        refinement_limits=refinement_limits,
        verbose=false, #  boolean in the output is produced or not
        depletion_handling=false,  # se true : motiplica epsilon_r per un fattore f nelle regioni non svuotate
        grid=Grid(sim05, max_tick_distance=max_tick_distance))

    println("💾 Saving simulation in $save_sim_path")
    ssd_write(save_sim_path, sim05)
end


final_plot = plot(
    plot(sim05.electric_potential, φ=0), # initial electric potential (boundary conditions)
    plot(sim05.point_types), # map of different point types: fixed point / inside or outside detector volume / depleted/undepleted
    plot(sim05.imp_scale), # impurity scale
    plot(sim05.q_eff_imp), # charge density distribution
    layout=(1, 4), size=(3000, 1000),
    title="grid tick $max_tick_distance and V_bias $V_bias "
)
savefig(final_plot, "plots/epot_problem/0.05mm_grid.png")


#--------------------- grid 0.02 ------------------

max_tick_distance = 0.02u"mm"
println("Eletric potential for grid $max_tick_distance")
# true  → ricalcola sempre e sovrascrivi
# false → usa il file salvato (se esiste)
recalculate = false   # <<<<<<<< CAMBIA QUI
save_sim_path = "saved_simulation/sim02"

if isfile(save_sim_path) && !recalculate
    println("⚡ Upoload simulation saved in : $save_sim_path")
    sim02 = ssd_read(save_sim_path, Simulation)

else
    println("🔧 New simulation for the electric potential...")

    sim02 = Simulation(SSD_examples[:InvertedCoax])
    sim02.detector = SolidStateDetector(sim02.detector, contact_id=2, contact_potential=V_bias)
    grid = Grid(sim02, max_tick_distance=max_tick_distance)


    # ATTENZIONE  A VEDERE SE IL POTENZIALE E' DEPLETED O NO
    calculate_electric_potential!(sim02,
        refinement_limits=refinement_limits,
        verbose=false, #  boolean in the output is produced or not
        depletion_handling=false,  # se true : motiplica epsilon_r per un fattore f nelle regioni non svuotate
        grid=Grid(sim02, max_tick_distance=max_tick_distance))

    println("💾 Saving simulation in $save_sim_path")
    ssd_write(save_sim_path, sim02)
end


final_plot = plot(
    plot(sim02.electric_potential, φ=0), # initial electric potential (boundary conditions)
    plot(sim02.point_types), # map of different point types: fixed point / inside or outside detector volume / depleted/undepleted
    plot(sim02.imp_scale), # impurity scale
    plot(sim02.q_eff_imp), # charge density distribution
    layout=(1, 4), size=(3000, 1000),
    title="grid tick $max_tick_distance and V_bias $V_bias "
)
savefig(final_plot, "plots/epot_problem/0.02mm_grid.png")

