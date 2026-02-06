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
using Printf

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
α = 2.0
fraction = [2.5, 2.0, 2.0]
P = prod(fraction)

# ===============================================
# Parametri di base
# ===============================================
#max_tick_fin = [0.05u"mm", 0.1u"rad", 0.5u"mm"]   # r, φ, z#
max_tick_fin = 0.5u"mm"            # oppure [0.01u"mm", 0.1u"rad", 0.01u"mm"]
max_tick_min = 5u"mm"               # oppure [5u"mm", 0u"rad", 5u"mm"]
α = 2.0
fraction = [2.5, 2.0, 2.0]
#fraction = [2.5, 2.0]   # fattori per il prodotto
P = prod(fraction)

# ===============================================
# Costruzione initial_tick coerente con max_tick_fin
# ===============================================
if max_tick_fin isa AbstractArray
    # 3 element array 
    println("max_tick_fin is an array")
    tmp = P .* max_tick_fin .* α
    initial_tick = [
        max(tmp[1], max_tick_min),
        tmp[2],  # φ: nessun minimo
        max(tmp[3], max_tick_min)
    ]
    initial_tick = Tuple(initial_tick)   # tuple r, φ, z
else
    # scalar
    println("max_tick_fin is a scalar")
    tmp = P .* max_tick_fin .* α

    initial_tick = P * max_tick_fin * α
    initial_tick = max(initial_tick, max_tick_min)
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

println("=== After the initial state, the max_tick_array for the refinement ===")
println("\nRefinement levels:")
for (i, ticks) in enumerate(max_tick_array)
    println("level $i → ", ticks)
end


# ============================================================
# INITIAL GRID + WP
# ============================================================
println("\n=== Starting Weighting Potential calculation with initial max tick = $initial_tick ===")
new_grid = Grid(sim, max_tick_distance=initial_tick)
apply_initial_state!(sim, WeightingPotential, 1, new_grid)
SolidStateDetectors.update_till_convergence!(sim, WeightingPotential, 1, depletion_handling=true, verbose=true)


println("Initial grid points: ",
    length(sim.weighting_potentials[1].grid.axes[1].ticks), " × ",
    length(sim.weighting_potentials[1].grid.axes[2].ticks), " × ",
    length(sim.weighting_potentials[1].grid.axes[3].ticks)
)



# ============================================================
# REFINEMENT DELLA GRIGLIA
# ============================================================
println("\n=== Starting grid refinement ===")
for ticks in max_tick_array
    # minima distanza sempre tuple di 3
    min_dist = (1e-5, 1e-5, 1e-5)

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
    @printf("max Δr = %.5g mm\n", maximum(diff(grid.axes[1].ticks)) * 1000)
    @printf("max Δz = %.5g mm\n", maximum(diff(grid.axes[3].ticks)) * 1000)
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

    println("r length: ", length(sim.weighting_potentials[1].grid.axes[1].ticks))
    println("ϕ length: ", length(sim.weighting_potentials[1].grid.axes[2].ticks))
    println("z length: ", length(sim.weighting_potentials[1].grid.axes[3].ticks))
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


if max_tick_fin isa AbstractArray
    r_mm = ustrip(u"mm", max_tick_fin[1])
    phi_r = ustrip(u"rad", max_tick_fin[2])
    z_mm = ustrip(u"mm", max_tick_fin[3])

    title_str = @sprintf(
        "max tick distance = r=%.3g mm, phi=%.3g rad, z=%.3g mm",
        r_mm, phi_r, z_mm
    )
else
    r_mm = ustrip(u"mm", max_tick_fin)
    phi_r = ustrip(u"rad", 0.1u"rad")
    z_mm = ustrip(u"mm", max_tick_fin)
    title_str = @sprintf(
        "max tick distance = r=%.3g mm, phi=%.3g rad, z=%.3g mm",
        r_mm, phi_r, z_mm
    )
end

x_min, x_max = 0.015, 0.020
y_min, y_max = 0.02, 0.03
x_rect = [x_min, x_max, x_max, x_min, x_min]
y_rect = [y_min, y_min, y_max, y_max, y_min]
p = plot(sim.weighting_potentials[1],
    contours_equal_potential=true,
    linecolor=:white,
    levels=5,
    title=title_str,
    size=(700, 400))
plot!(sim.detector, st=:slice, φ=90u"°", legend=false)
plot!(x_rect, y_rect,
    seriestype=:shape,
    linecolor=:white,
    lw=1.5,
    fillalpha=0,
    label=""
)

fname = @sprintf(
    "ref_max_tick_r%.3gmm__z%.3gmm.png",
    r_mm, z_mm
)

savefig(p, joinpath(
    "notebooks/new_refinement_on_grid/vector_3_refinement",
    fname
))
