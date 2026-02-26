cd("/home/ritaferi/Phd/SSDdev/")
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using SolidStateDetectors
using Unitful
using Plots
using LegendHDF5IO
using HDF5

gr()

refinement_limits = [0.2, 0.1, 0.05, 0.02]


save_sim_path = "saved_simulation/sim_ILM.h5"


recalculate = true   # <<<<<<<< CAMBIA QUI
if isfile(save_sim_path) && !recalculate
    println("⚡ Upoload simulation saved in : $save_sim_path")
    sim = ssd_read(save_sim_path, Simulation)

else
    println("🔧 New simulation for the electric potential...")
    max_tick_distance_ep = 0.08u"mm"
    #sim = Simulation(SSD_examples[:InvertedCoax])
    sim = Simulation(SSD_examples[:IVCIlayer])
    sim.detector = SolidStateDetector(sim.detector, contact_id=2, contact_potential=500u"V")
    grid = Grid(sim)


    calculate_electric_potential!(sim,
        refinement_limits=[0.2, 0.1, 0.05, 0.02],
        verbose=true, #  boolean in the output is produced or not
        depletion_handling=true,  # motiplica epsilon_r per un fattore f nelle regioni non svuotate
        grid=grid)

    println("💾 Saving simulation in $save_sim_path")
    ssd_write(save_sim_path, sim)
end



refinement_limits = [0.2, 0.1, 0.05, 0.02]
#=
println("======== Electric Potential calculation =======")

sim = Simulation(SSD_examples[:IVCIlayer])
sim.detector = SolidStateDetector(sim.detector, contact_id=2, contact_potential=500u"V")
SolidStateDetectors._calculate_potential_max_tick_refinement!(sim,
    ElectricPotential,
    depletion_handling=true,
    verbose=true,
    refinement_limits=refinement_limits)

pel = plot(sim.electric_potential,
    contours_equal_potential=true, full_det=true, linecolor=:white,
    levels=5,
    size=(700, 400),
    legend=false)
plot!(sim.detector, full_det=true, st=:slice, φ=0)

savefig(pel, "electric_pot_new_fnc.png")

#max_tick_distance = (0.5u"mm", 1u"rad", 0.5u"mm")
#max_tick_distance = 0.05u"mm"

grid_ep = sim.electric_potential.grid
println("max Δr =$(maximum(diff(grid_ep.axes[1].ticks)) * 1000) mm")
println("max Δz = $(maximum(diff(grid_ep.axes[3].ticks)) * 1000) mm")

=#


println("======== Weighting Potential calculation =======")

#max_tick_distance = (4u"mm", 1u"rad", 0.6u"mm")

max_tick_distance = 1u"mm"

ref_grid = [2.5]
# inizio del time
t0 = time_ns()

# passo la griglia del pot_elettrico
SolidStateDetectors._calculate_potential_max_tick_refinement!(sim,
    WeightingPotential,
    1,
    max_tick_distance=max_tick_distance,
    grid=sim.electric_potential.grid,
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

#key_name = "grid_$(ustrip(max_tick_distance))mm"

#=
# --- dati da salvare ---
time_seconds = elapsed
min_potential = min_wp



# directory di salvataggio
save_dir = "saved_results_new_function"
mkpath(save_dir)   # crea la directory se non esiste

# nome file
filename = joinpath(save_dir, "simulation_with_$(1)_jumps.h5")


file_mode = isfile(filename) ? "r+" : "w"
# apertura in append per non sovrascrivere le altre simulazioni
h5open(filename, file_mode) do f
    # se il gruppo esiste già, lo riutilizziamo, altrimenti lo creiamo
    g = haskey(f, key_name) ? f[key_name] : create_group(f, key_name)

    # scrivi i dataset
    g["time_s"] = time_seconds
    g["min_weighting_potential"] = min_potential
end

println("💾 File salvato in: $(abspath(filename))")
println("   Gruppo: $key_name")
 =#

#save_dir_plot = "saved_plots_new_function_only_initialized"
#mkpath(save_dir_plot)
x_min, x_max = 0.015, 0.020
y_min, y_max = 0.02, 0.03
x_rect = [x_min, x_max, x_max, x_min, x_min]
y_rect = [y_min, y_min, y_max, y_max, y_min]
p = plot(sim.weighting_potentials[1],
    contours_equal_potential=true,
    linecolor=:white,
    levels=5,
    title="max tick Ep= $max_tick_distance_ep",
    size=(700, 400))
plot!(sim.detector, st=:slice, φ=90u"°", legend=false)
plot!(x_rect, y_rect,
    seriestype=:shape,
    linecolor=:white,
    lw=1.5,
    fillalpha=0,
    label=""
)
#fname = "only_1ref_$max_tick_distance.png"
#=
savefig(p, joinpath(
    save_dir_plot,
    fname
))=#

#savefig(p, "testing_grid_from_Ep.png")
