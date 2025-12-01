# IN QUESTO CODICE VALUTO IL COMPORTAMENTO DEL POTENZIALE PESATO AL VARIARE 
# DEI DIVERSI PARAMETRI DI INTERESSE

cd("/home/ritaferi/Phd/SSDdev")
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using SolidStateDetectors
using Unitful
using Plots
using LegendHDF5IO

gr()

max_tick_distance = 0.1u"mm"
refinement_limits = [0.2, 0.1, 0.05, 0.02]
save_sim_path = "saved_simulation/sim.h5"



"""
calculate_electric_potential!(sim,
    refinement_limits=refinement_limits,
    verbose=false, #  boolean in the output is produced or not
    depletion_handling=true,  # motiplica epsilon_r per un fattore f nelle regioni non svuotate
    grid=grid)
ssd_write("saved_simulation/sim.h5", sim)
# println("initial grid is the same of the final grid for electric potential: ", isequal(grid, grid_for_Ep) ? "Yes" : "No")

sim = ssd_read("notebooks/sim.h5", Simulation)
"""


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

    sim = Simulation(SSD_examples[:InvertedCoax])
    sim.detector = SolidStateDetector(sim.detector, contact_id=2, contact_potential=500u"V")
    grid = Grid(sim, max_tick_distance=max_tick_distance)


    calculate_electric_potential!(sim,
        refinement_limits=refinement_limits,
        verbose=false, #  boolean in the output is produced or not
        depletion_handling=true,  # motiplica epsilon_r per un fattore f nelle regioni non svuotate
        grid=grid)

    println("💾 Saving simulation in $save_sim_path")
    ssd_write(save_sim_path, sim)
end

println("Starting weighting potential simulation ....")


calculate_weighting_potential!(sim, 1,
    max_n_iterations=-1,
    convergence_limit=5e-7,   # PARAMETRO DA CAMBIARE
    refinement_limits=refinement_limits,
    depletion_handling=true,
    max_tick_distance=max_tick_distance,
    grid=grid = Grid(sim, for_weighting_potential=true,
        max_tick_distance=max_tick_distance))







"""
for contact in sim.detector.contacts
    calculate_weighting_potential!(sim, contact.id,
        refinement_limits=[0.2, 0.1, 0.05, 0.02],
        verbose=false,
        depletion_handling=true,
        grid=Grid(sim, max_tick_distance=max_tick_distance, for_weighting_potential=true)
    )
end



grd_el1 = sim.weighting_potentials[1].grid
grd_el2 = sim.weighting_potentials[2].grid


println("compare if electric and weighting potential have the same grid: ", isequal(grid_for_Ep, grd_el1) ? "Yes" : "No")


wp1 = plot(sim.weighting_potentials[1], contours_equal_potential=true,
    linecolor=:white, levels=5)
plot!(sim.detector, st=:slice, φ=0)
wp2 = plot(sim.weighting_potentials[2], contours_equal_potential=true,
    linecolor=:white, levels=5)
plot!(sim.detector, st=:slice, φ=0)
w = plot(
    wp1, wp2,
    size=(900, 700), layout=(1, 2)
)
"""

wp1 = plot(sim.weighting_potentials[1], contours_equal_potential=true,
    linecolor=:white, levels=5)
plot!(sim.detector, st=:slice, φ=0, legend=false)

#  calcolo il potenziale solo per il 1 elettrodo (quello problematico)
savefig(wp1, "plots/0.1mm_5e-7_convergence_limit.png")
