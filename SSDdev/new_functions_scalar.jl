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

#------------ qui inizia la simulazione del weighting potential -------
# inizializzo lo stato con una griglia più grande partendo dai default parameter
# probabilmente partire dal max tick di defautl e poi raffinare suil max_tick è troppo. ( passiamo da 20 mm a 0.5 mm)
# provo a partire da 2 il primo max tick e poi raffinare

max_tick_fin = 0.03u"mm" # oppure max_tick_fin = [0.01u"mm", 1u"rad", 0.01u"mm"]
max_tick_min = 5u"mm"  # oppure max_tick_min = [51u"mm", 1u"rad", 5u"mm"]

α = 2.0
fraction = [2.5, 2, 2]

P = prod(fraction)
initial_tick = P * max_tick_fin * α
initial_tick = max(initial_tick, max_tick_min)


println("grid_final = ", max_tick_fin)
println("Prodotto refinement = ", P)
println("grid_init = ", initial_tick)

max_tick_array = (max_tick_fin * fraction[1] * fraction[2] * fraction[3], max_tick_fin * fraction[1] * fraction[2],
    max_tick_fin * fraction[1],
    max_tick_fin)

println("max tick array = ", max_tick_array)


#max_tick_array = (0.5u"mm", 0.35u"mm", 0.2u"mm", 0.1u"mm")
println("=== Starting Weighting Potential calculation with initial max tick = $initial_tick ===")
new_grid = Grid(sim, max_tick_distance=initial_tick)
apply_initial_state!(sim, WeightingPotential, 1, new_grid) # optional
SolidStateDetectors.update_till_convergence!(sim, WeightingPotential, 1, depletion_handling=true, verbose=true)
println("Initial grid initial state: ", length(sim.weighting_potentials[1].grid.axes[1].ticks))



"""
inizio con il refinement della griglia 
"""
refinement_limits = Float32.([0.2, 0.1, 0.05, 0.02])
n_refinement_steps = length(refinement_limits)



for max_tick in max_tick_array
    println(" Applying max tick distance = $max_tick ")
    SolidStateDetectors.refine_max_tick!(sim, WeightingPotential, 1, (max_tick, 1u"rad", max_tick),
        (1e-15u"mm", 1e-15u"rad", 1e-15u"mm")) # qui metterò i valori con la fuzione di claudia
    SolidStateDetectors.update_till_convergence!(sim, WeightingPotential, 1, depletion_handling=true, verbose=true)
    println("r length: ", length(sim.weighting_potentials[1].grid.axes[1].ticks))
    println("ϕ length: ", length(sim.weighting_potentials[1].grid.axes[2].ticks))
    println("z length: ", length(sim.weighting_potentials[1].grid.axes[3].ticks))

    grid = sim.weighting_potentials[1].grid
    println("max Δr = ", maximum(diff(grid.axes[1].ticks)), " m")
    println("max Δz = ", maximum(diff(grid.axes[3].ticks)), " m")

end



println("================== Starting refinement on the potential ================")





# faccio il refinement sul potenziale solo per l'ultimo max_ticxk
for iref in 1:n_refinement_steps
    println(" Refinement step $(refinement_limits[iref])")
    max_diff_array = (refinement_limits[iref], refinement_limits[iref], refinement_limits[iref])
    SolidStateDetectors.refine!(sim, WeightingPotential, 1, max_diff_array)
    SolidStateDetectors.update_till_convergence!(sim, WeightingPotential, 1, depletion_handling=true, verbose=true)
    println("r length: ", length(sim.weighting_potentials[1].grid.axes[1].ticks))
    println("ϕ length: ", length(sim.weighting_potentials[1].grid.axes[2].ticks))
    println("z length: ", length(sim.weighting_potentials[1].grid.axes[3].ticks))
end




wp = sim.weighting_potentials[1]

r_ticks = wp.grid.axes[1].ticks   # r
z_ticks = wp.grid.axes[3].ticks   # z

rmin, rmax = ustrip.(u"m", (0.015u"m", 0.02u"m"))
zmin, zmax = ustrip.(u"m", (0.02u"m", 0.03u"m"))

r_inds = findall(r -> rmin ≤ r ≤ rmax, r_ticks)
z_inds = findall(z -> zmin ≤ z ≤ zmax, z_ticks)

isempty(r_inds) && error("No points found in the r interval")
isempty(z_inds) && error("No points found in the z interval")

vals = wp.data[r_inds, 1, z_inds]
min_wp = minimum(vals)

println(">>> Min WeightingPotential in the test volume = $min_wp")






x_min, x_max = 0.015, 0.020
y_min, y_max = 0.02, 0.03
x_rect = [x_min, x_max, x_max, x_min, x_min]
y_rect = [y_min, y_min, y_max, y_max, y_min]
p = plot(sim.weighting_potentials[1],
    contours_equal_potential=true,
    linecolor=:white,
    levels=5,
    title="max tick distance = $(max_tick_array[end])",
    size=(400, 400))
plot!(sim.detector, st=:slice, φ=90u"°", legend=false)
plot!(x_rect, y_rect,
    seriestype=:shape,
    linecolor=:white,
    lw=1.5,
    fillalpha=0,
    label=""
)

savefig(p, "notebooks/new_refinement_on_grid/ne_plots_correceted/ref_max_tick_$max_tick_fin small initial grid.png")

