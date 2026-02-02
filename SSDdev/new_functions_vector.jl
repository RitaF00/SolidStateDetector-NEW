# ============================================================
# SETUP AMBIENTE
# ============================================================
cd("/home/ritaferi/Phd/SSDdev/")

using Pkg
Pkg.activate(".")
Pkg.instantiate()

using SolidStateDetectors
using Unitful
using Plots
using LegendHDF5IO

gr()  # backend per Plots

# ============================================================
# PARAMETRI
# ============================================================
refinement_limits = Float32.([0.2, 0.1, 0.05, 0.02])
save_sim_path = "saved_simulation/sim_ILM.h5"
recalculate = false   # <<< CAMBIA QUI

# ============================================================
# CARICAMENTO / CALCOLO DELLA SIMULAZIONE ELETTRICA
# ============================================================
if isfile(save_sim_path) && !recalculate
    println("⚡ Upload simulation saved in : $save_sim_path")
    sim = ssd_read(save_sim_path, Simulation)
else
    println("🔧 New simulation for the electric potential...")
    max_tick_distance = 0.5u"mm"
    sim = Simulation(SSD_examples[:IVCIlayer])
    sim.detector = SolidStateDetector(sim.detector, contact_id=2, contact_potential=500u"V")
    grid = Grid(sim, max_tick_distance=max_tick_distance)

    calculate_electric_potential!(sim,
        refinement_limits=refinement_limits,
        verbose=true,
        depletion_handling=true,
        grid=grid
    )

    println("💾 Saving simulation in $save_sim_path")
    ssd_write(save_sim_path, sim)
end

# ============================================================
# SETUP GRID PER WEIGHTING POTENTIAL
# ============================================================
# max_tick scalare o array
max_tick_fin = [0.01u"mm", 0.1u"rad", 0.01u"mm"]   # r, φ, z
max_tick_fin = 0.01u"mm"


α = 2.0
fraction = [2.5, 2.0, 2.0]
P = prod(fraction)

# ===============================================
# Parametri di base
# ===============================================
max_tick_fin = 0.01u"mm"            # oppure [0.01u"mm", 0.1u"rad", 0.01u"mm"]
max_tick_min = 5u"mm"               # oppure [5u"mm", 0u"rad", 5u"mm"]
α = 2.0
fraction = [2.5, 2.0, 2.0]     # fattori per il prodotto
P = prod(fraction)

# ===============================================
# Costruzione initial_tick coerente con max_tick_fin
# ===============================================
if max_tick_fin isa AbstractArray
    # array di 3 elementi → initial_tick come array di 3 elementi
    initial_tick = max.(P .* max_tick_fin .* α, max_tick_min)
    initial_tick_for_grid = Tuple(initial_tick)   # tuple r, φ, z
else
    # scalare → initial_tick scalare
    initial_tick = max(P * max_tick_fin * α, max_tick_min)
    initial_tick_for_grid = (initial_tick, 0.1u"rad", initial_tick)         # rimane scalare
end

println("grid_final = ", max_tick_fin)
println("Prodotto refinement = ", P)
println("grid_init = ", initial_tick)

# ============================================================
# COSTRUZIONE MAX TICK ARRAY (livelli di refinement) – GENERALIZZATA
# ============================================================

levels = reverse(vcat(1.0, cumprod(fraction)))  # fattori per i livelli

# funzione interna per normalizzare max_tick in tuple di 3 o scalar
normalize_max_tick(x) = begin
    if x isa Number || x isa Quantity
        # scalar → ritorna Quantity scalare
        return x
    elseif length(x) == 3
        # tuple o array di 3 → ritorna tuple di 3
        return Tuple(x)
    else
        error("max_tick deve essere scalar o una collezione di 3 elementi")
    end
end

# costruzione max_tick_array generalizzata
max_tick_array = [normalize_max_tick(level .* max_tick_fin) for level in levels]

println("\nRefinement levels:")
for (i, ticks) in enumerate(max_tick_array)
    println("level $i → ", ticks)
end


# ============================================================
# INITIAL GRID + WP
# ============================================================
println("\n=== Starting Weighting Potential calculation with initial max tick = $initial_tick ===")
new_grid = Grid(sim, max_tick_distance=initial_tick_for_grid)
apply_initial_state!(sim, WeightingPotential, 1, new_grid)
SolidStateDetectors.update_till_convergence!(sim, WeightingPotential, 1, depletion_handling=true, verbose=true)

grid = sim.weighting_potentials[1].grid
println("Initial grid points: ",
    length(grid.axes[1].ticks), " × ",
    length(grid.axes[2].ticks), " × ",
    length(grid.axes[3].ticks)
)

# plot iniziale
p = plot(sim.weighting_potentials[1], size=(400, 400))
plot!(sim.detector, st=:slice, φ=90u"°", legend=false)
savefig(p, "initial_WP_state.png")

# ============================================================
# REFINEMENT DELLA GRIGLIA
# ============================================================
println("\n=== Starting grid refinement ===")
for ticks in max_tick_array
    # minima distanza sempre tuple di 3
    min_dist = (1e-15u"mm", 1e-15u"rad", 1e-15u"mm")

    println("\nApplying max tick distance = ", ticks)

    # se ticks è scalar, rimane scalar, se è tuple di 3, rimane tuple
    SolidStateDetectors.refine_max_tick!(sim, WeightingPotential, 1, ticks, min_dist)

    SolidStateDetectors.update_till_convergence!(
        sim, WeightingPotential, 1,
        depletion_handling=true,
        verbose=true
    )

    grid = sim.weighting_potentials[1].grid
    println("r length: ", length(grid.axes[1].ticks))
    println("z length: ", length(grid.axes[3].ticks))
    println("max Δr = ", maximum(diff(grid.axes[1].ticks)))
    println("max Δz = ", maximum(diff(grid.axes[3].ticks)))
end


# ============================================================
# REFINEMENT SUL POTENZIALE
# ============================================================
println("\n================== Starting refinement on the potential ==================")
n_refinement_steps = length(refinement_limits)
for iref in 1:n_refinement_steps
    println("\nRefinement step with limit = ", refinement_limits[iref])
    max_diff_array = max_diff_array = (refinement_limits[iref], refinement_limits[iref], refinement_limits[iref])
    SolidStateDetectors.refine!(sim, WeightingPotential, 1, max_diff_array)
    SolidStateDetectors.update_till_convergence!(sim, WeightingPotential, 1, depletion_handling=true, verbose=true)

    grid = sim.weighting_potentials[1].grid
    println("r length: ", length(grid.axes[1].ticks))
    println("ϕ length: ", length(grid.axes[2].ticks))
    println("z length: ", length(grid.axes[3].ticks))
end

# ============================================================
# VERIFICA VALORI DEL WP NEL VOLUME DI TEST
# ============================================================
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

# ============================================================
# PLOT CON RETTANGOLO DI CONTROLLO
# ============================================================
x_min, x_max = 0.015, 0.020
y_min, y_max = 0.02, 0.03
x_rect = [x_min, x_max, x_max, x_min, x_min]
y_rect = [y_min, y_min, y_max, y_max, y_min]

p = plot(
    sim.weighting_potentials[1],
    contours_equal_potential=true,
    linecolor=:white,
    levels=5,
    title="max tick distance = $(max_tick_array[end, :])",
    size=(400, 400)
)
plot!(sim.detector, st=:slice, φ=90u"°", legend=false)
plot!(x_rect, y_rect,
    seriestype=:shape,
    linecolor=:white,
    lw=1.5,
    fillalpha=0,
    label=""
)

savefig(p, "notebooks/new_refinement_on_grid/plot_new_refinement/ref_max_tick.png")
