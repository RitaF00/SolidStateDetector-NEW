cd("/home/ritaferi/Phd/SSDdev/")

using Pkg
Pkg.activate(".")
Pkg.instantiate()

using SolidStateDetectors
using Unitful
using Plots
using LegendHDF5IO
using Printf


gr()
refinement_limits = Float32.([0.2, 0.1, 0.05, 0.02])
save_sim_path = "saved_simulation/sim_ILM.h5"
recalculate = true

if isfile(save_sim_path) && !recalculate
    println("⚡ Upload simulation saved in : $save_sim_path")
    sim = ssd_read(save_sim_path, Simulation)
else
    println("🔧 New simulation for the electric potential...")
    max_tick_distance = 0.5u"mm"
    sim = Simulation(SSD_examples[:IVCIlayer])
    sim.detector = SolidStateDetector(sim.detector, contact_id=2, contact_potential=500u"V")
    calculate_electric_potential!(sim,
        refinement_limits=refinement_limits,
        verbose=true,
        depletion_handling=true,
        #max_tick_distance=max_tick_distance
    )

    println("💾 Saving simulation in $save_sim_path")
    ssd_write(save_sim_path, sim)
end

α = 2.0
fraction = [2.5, 2.0, 2.0]
grid_ep = sim.electric_potential.grid
println("Electric potential grid points: ",
    length(grid_ep.axes[1].ticks), " × ",
    length(grid_ep.axes[2].ticks), " × ",
    length(grid_ep.axes[3].ticks)
)



println("=============== Applying initial state =======================")
apply_initial_state!(sim, WeightingPotential, 1, grid_ep)
SolidStateDetectors.update_till_convergence!(sim, WeightingPotential, 1, depletion_handling=true, verbose=true)
println("Initial grid points: ",
    length(sim.weighting_potentials[1].grid.axes[1].ticks), " × ",
    length(sim.weighting_potentials[1].grid.axes[2].ticks), " × ",
    length(sim.weighting_potentials[1].grid.axes[3].ticks)
)



println("\n================== Starting refinement on the potential ==================")
n_refinement_steps = length(refinement_limits)
for iref in 1:n_refinement_steps
    println("\nRefinement step with limit = ", refinement_limits[iref])
    max_diff_array = (refinement_limits[iref], refinement_limits[iref], refinement_limits[iref])
    SolidStateDetectors.refine!(sim, WeightingPotential, 1, max_diff_array)
    SolidStateDetectors.update_till_convergence!(sim, WeightingPotential, 1, depletion_handling=true, verbose=true)

    println("r length: ", length(sim.weighting_potentials[1].grid.axes[1].ticks))
    println("ϕ length: ", length(sim.weighting_potentials[1].grid.axes[2].ticks))
    println("z length: ", length(sim.weighting_potentials[1].grid.axes[3].ticks))
end

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

x_min, x_max = 0.015, 0.020
y_min, y_max = 0.02, 0.03
x_rect = [x_min, x_max, x_max, x_min, x_min]
y_rect = [y_min, y_min, y_max, y_max, y_min]
p = plot(sim.weighting_potentials[1],
    contours_equal_potential=true,
    linecolor=:white,
    levels=5,
    title="ep initial grid ,
    $(length(sim.weighting_potentials[1].grid.axes[1].ticks)) ×  $(length(sim.weighting_potentials[1].grid.axes[2].ticks)) × $(length(sim.weighting_potentials[1].grid.axes[3].ticks))",
    size=(1200, 900))
plot!(sim.detector, st=:slice, φ=90u"°", legend=false)
plot!(x_rect, y_rect,
    seriestype=:shape,
    linecolor=:white,
    lw=1.5,
    fillalpha=0,
    label=""
)

savefig(
    p,
    "plot_custom_grid/electric_pot_grig$(length(sim.weighting_potentials[1].grid.axes[1].ticks)), × ,$(length(sim.weighting_potentials[1].grid.axes[2].ticks)), × , $(length(sim.weighting_potentials[1].grid.axes[3].ticks)).png"
)


println("ok")