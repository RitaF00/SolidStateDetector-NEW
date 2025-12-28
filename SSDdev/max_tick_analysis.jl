cd("/home/ritaferi/Phd/SSDdev/")
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using SolidStateDetectors
using Unitful
using Plots
using LegendHDF5IO

gr()



#------------ definizione dei parametri ------

max_tick_distance = 0.5u"mm"
refinement_limits = [0.2, 0.1, 0.05, 0.02]

#-------- ICPC vuoto -----------
save_sim_path = "saved_simulation/sim.h5"
#-------- ICPC ILM -----------
#save_sim_path = "saved_simulation/sim_ILM.h5"



# ------------------------
# SCELTA 
# ------------------------
# true  → ricalcola sempre e sovrascrivi
# false → usa il file salvato (se esiste)
recalculate = false   # <<<<<<<< CAMBIA QUI

# ------------------------
# Logica principale
# ------------------------
if isfile(save_sim_path) && !recalculate
    println("⚡ Upoload simulation saved in : $save_sim_path")
    sim = ssd_read(save_sim_path, Simulation)

else
    println("🔧 New simulation for the electric potential...")

    #sim = Simulation(SSD_examples[:InvertedCoax])
    sim = Simulation(SSD_examples[:IVCIlayer])
    sim.detector = SolidStateDetector(sim.detector, contact_id=2, contact_potential=500u"V")
    grid = Grid(sim, max_tick_distance=max_tick_distance)


    calculate_electric_potential!(sim,
        refinement_limits=[0.2, 0.1, 0.05, 0.02],
        verbose=false, #  boolean in the output is produced or not
        depletion_handling=true,  # motiplica epsilon_r per un fattore f nelle regioni non svuotate
        grid=grid)

    println("💾 Saving simulation in $save_sim_path")
    ssd_write(save_sim_path, sim)
end

println("Simulating weighting potential ")


max_tick_distance_array = [0.5u"mm", 0.45u"mm", 0.4u"mm", 0.35u"mm",
    0.3u"mm", 0.25u"mm", 0.2u"mm", 0.15u"mm",
    0.1u"mm"]



n_rows, n_cols = 2, 5
plot_list = []

#for max_tick_distance in max_tick_distance_array
#println("Simulating weighting potential with max tickness = $max_tick_distance")

# Calcolo del weighting potential solo per il primo elettrodo
calculate_weighting_potential!(sim, 1,
    refinement_limits=refinement_limits,
    depletion_handling=true,
    grid=Grid(sim,
        for_weighting_potential=true))#,
#max_tick_distance=0.3u"mm")) # sto dacendo default


# Creiamo il plot
p = plot(sim.weighting_potentials[1],
    contours_equal_potential=true,
    linecolor=:white, levels=5)
#title="max_tick = $(max_tick_distance)")
plot!(sim.detector, st=:slice, φ=0, legend=false)

#    push!(plot_list, p)
#end


#final_plot = plot(plot_list..., layout=(n_rows, n_cols), size=(2000, 800))

savefig(p, "plots/max_tick-analysis/default_vacuum.png")


