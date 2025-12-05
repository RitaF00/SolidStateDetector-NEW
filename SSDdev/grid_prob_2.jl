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
#refinement_limits = [0.2, 0.1, 0.05, 0.02]
refinement_limits = [0.2, 0.1, 0.05]
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


# Lista dei convergence limit da testare
convergence_limits = [1e-7]

# Creiamo una figura vuota con 2 righe e 5 colonne
n_rows, n_cols = 1, 1
plot_list = []

for cl in convergence_limits
    println("Simulating weighting potential with convergence_limit = $cl")

    # Generiamo una nuova grid per ogni calcolo (puoi eventualmente riusare la stessa)
    grid_wp = Grid(sim, for_weighting_potential=true, max_tick_distance=max_tick_distance)

    # Calcolo del weighting potential solo per il primo elettrodo
    calculate_weighting_potential!(sim, 1,
        max_n_iterations=-1,  # non faccio continuaare all'infinito
        convergence_limit=cl,
        refinement_limits=refinement_limits,
        depletion_handling=true,
        max_tick_distance=max_tick_distance,
        n_iterations_between_checks=30000,
        grid=grid_wp)

    # Creiamo il plot
    p = plot(sim.weighting_potentials[1],
        contours_equal_potential=true,
        linecolor=:white, levels=5,
        title="cl = $(cl)")
    plot!(sim.detector, st=:slice, φ=0, legend=false)

    push!(plot_list, p)
end

# Creiamo un'unica figura con layout 2x5
final_plot = plot(plot_list..., layout=(n_rows, n_cols), size=(2000, 800))
savefig(final_plot, "plots/update_convergence_plot/30_000_bewtween_checks_1e-7_no_0.02_ref.png")







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


wp1 = plot(sim.weighting_potentials[1], contours_equal_potential=true,
    linecolor=:white, levels=5)
plot!(sim.detector, st=:slice, φ=0, legend=false)

#  calcolo il potenziale solo per il 1 elettrodo (quello problematico)
savefig(wp1, "plots/convergences/0.1mm_5e-8_convergence_limit.png")
"""