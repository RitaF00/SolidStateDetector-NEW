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
    max_tick_distance = 0.5u"mm"
    sim = Simulation(SSD_examples[:InvertedCoax])
    #sim = Simulation(SSD_examples[:IVCIlayer])
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


Δz = 0.25u"mm"
f_conversion = sim.world.intervals[3].right / sim.world.intervals[1].right
# devi ricordarti che dal file yaml le lunghezze vengono espresse in metri
Δr = sim.world.intervals[1].right * 1000u"mm" / 4
#max_tick_distance = (Δr, 0u"rad", Δz)


#max_tick_distance = 0.25u"mm"
n_rows, n_cols = 2, 5
plot_list = []

max_tick_array = [0.5u"mm", 0.45u"mm", 0.4u"mm", 0.35u"mm", 0.3u"mm", 0.25u"mm", 0.2u"mm", 0.15u"mm", 0.1u"mm"]

#=
# Calcolo del weighting potential solo per il primo elettrodo
for max_tick_distance in max_tick_array
    println(" max_tick_distance = $max_tick_distance ")
    calculate_weighting_potential!(sim, 1,
        refinement_limits=refinement_limits,
        depletion_handling=true,
        grid=Grid(sim,
            for_weighting_potential=true,
            max_tick_distance=max_tick_distance))
end

=#


max_tick_distance = 0.5u"mm"

max_tick_array = [0.5u"mm", 0.45u"mm", 0.4u"mm", 0.35u"mm", 0.3u"mm", 0.25u"mm", 0.2u"mm", 0.15u"mm", 0.1u"mm", 0.05u"mm", 0.02u"mm"]



# Assumiamo che max_tick_array esista già
n = length(max_tick_array)

# Creiamo una griglia di subplot (ad esempio 2 colonne)
ncols = 3
nrows = ceil(Int, n / ncols)

# Preparo il layout
plt = plot(layout=(nrows, ncols), size=(800, 400 * nrows))

# Coordinate del rettangolo
x_min = 0.015
x_max = 0.020
y_min = 0.02
y_max = 0.04
x_rect = [x_min, x_max, x_max, x_min, x_min]
y_rect = [y_min, y_min, y_max, y_max, y_min]

# Loop su tutti i weighting potentials
for (i, max_tick_distance) in enumerate(max_tick_array)

    println(" calculating weighting potential with max_tick_distance = $max_tick_distance ")

    # Se vuoi, puoi calcolarlo qui o usare quelli già calcolati
    calculate_weighting_potential!(sim, 1,
        refinement_limits=refinement_limits,
        depletion_handling=true,
        grid=Grid(sim,
            for_weighting_potential=true,
            max_tick_distance=max_tick_distance))

    # Plot principale
    plot!(plt[i], sim.weighting_potentials[1],
        contours_equal_potential=true,
        linecolor=:white,
        levels=5,
        title="max tick distance = $max_tick_distance")

    # Plot del detector
    plot!(plt[i], sim.detector, st=:slice, φ=0, legend=false)

    # Rettangolo
    plot!(plt[i], x_rect, y_rect,
        seriestype=:shape,
        linecolor=:white,
        lw=1.5,
        fillalpha=0,
        label="")
end

display(plt)
savefig(plt, "weighting_potential_vs_max_tick_distance.png")


