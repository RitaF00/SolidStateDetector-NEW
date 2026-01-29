cd("/home/ritaferi/Phd/SSDdev/")
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using SolidStateDetectors
using Unitful
using Plots
using LegendHDF5IO

gr()




#-------- ICPC vuoto -----------
#save_sim_path = "saved_simulation/sim.h5"
#-------- ICPC ILM -----------
save_sim_path = "saved_simulation/sim_ILM.h5"



# ------------------------
# SCELTA 
# ------------------------
# true  → ricalcola sempre e sovrascrivi
# false → usa il file salvato (se esiste)
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

    refinement_limits = [0.2, 0.1, 0.05, 0.02]
    calculate_electric_potential!(sim,
        refinement_limits=[0.2, 0.1, 0.05, 0.02],
        verbose=true, #  boolean in the output is produced or not
        depletion_handling=true,  # motiplica epsilon_r per un fattore f nelle regioni non svuotate
        grid=grid)

    println("💾 Saving simulation in $save_sim_path")
    ssd_write(save_sim_path, sim)
end

#------------ qui inizia la simulazione del weighting potential -------
# inizializzo lo stato con una griglia più grande con i parametri di default


#=

p = plot(sim.weighting_potentials[1], size=(400, 400))
plot!(sim.detector, st=:slice, φ=90u"°", legend=false)
savefig(p, "initial_WP_state.png")



"""
inizio con il refinement della griglia costruita con i defaul parameters
"""
refinement_limits = Float32.([0.2, 0.1, 0.05, 0.02])
n_refinement_steps = length(refinement_limits)
for iref in 1:n_refinement_steps
    println("Refinement step $(refinement_limits[iref])")
    max_diff_array = (refinement_limits[iref], refinement_limits[iref], refinement_limits[iref])
    SolidStateDetectors.refine!(sim, WeightingPotential, 1, max_diff_array)
    SolidStateDetectors.update_till_convergence!(sim, WeightingPotential, 1, depletion_handling=true, verbose=true)
end

d = 0.5u"mm"

#=
sim.weighting_potentials[1] =
    SolidStateDetectors.refine_scalar_potential_by_tick_distance(
        sim.weighting_potentials[1],
        (d, d, d),
        (0.0001u"mm", 0.0001u"mm", 0.0001u"mm")
    );
=#

"""
appplico il refine max tick direttamente sull'array
"""

SolidStateDetectors.refine_max_tick!(sim, WeightingPotential, 1, (d, d, d),
    (1e-8u"mm", 1e-8u"mm", 1e-8u"mm"))
SolidStateDetectors.update_till_convergence!(sim, WeightingPotential, 1, depletion_handling=true, verbose=true)

"""
applico nuovamente il refinement
"""
refinement_limits = Float32.([0.2, 0.1, 0.5, 0.2])
n_refinement_steps = length(refinement_limits)
for iref in 1:n_refinement_steps
    println("Refinement step $(refinement_limits[iref])")
    max_diff_array = (refinement_limits[iref], refinement_limits[iref], refinement_limits[iref])
    SolidStateDetectors.refine!(sim, WeightingPotential, 1, max_diff_array)
    SolidStateDetectors.update_till_convergence!(sim, WeightingPotential, 1, depletion_handling=true, verbose=true)
end
=#


#max_tick_array = [0.5u"mm", 0.4u"mm", 0.3u"mm", 0.2u"mm", 0.1u"mm"]


fraction = [2.5, 2, 2]


max_tick_fin = 0.1u"mm"

# applico la orima griglia con i max tick di default
new_grid = Grid(sim)
apply_initial_state!(sim, WeightingPotential, 1, new_grid) # optional
SolidStateDetectors.update_till_convergence!(sim, WeightingPotential, 1, depletion_handling=true, verbose=true)
println("Initial grid initial state: ", length(sim.weighting_potentials[1].grid.axes[1].ticks))

println("===================================================")
println("max_tick_fin = $max_tick_fin")

# Costruzione gerarchica dei max_tick
max_tick_array = (
    max_tick_fin * fraction[1] * fraction[2] * fraction[3],
    max_tick_fin * fraction[1] * fraction[2],
    max_tick_fin * fraction[1],
    max_tick_fin
)

for max_tick in max_tick_array
    println("  Applying max tick distance = $max_tick")
    # qui faccio il refinemt sulla griglia
    SolidStateDetectors.refine_max_tick!(
        sim,
        WeightingPotential,
        1,
        (max_tick, max_tick, max_tick),
        (1e-15u"mm", 1e-15u"mm", 1e-15u"mm")
    )

    SolidStateDetectors.update_till_convergence!(
        sim,
        WeightingPotential,
        1,
        depletion_handling=true,
        verbose=true
    )

    refinement_limits = Float32.([0.2, 0.1, 0.05, 0.02])
    n_refinement_steps = length(refinement_limits)
    for iref in 1:n_refinement_steps
        println("Refinement step $(refinement_limits[iref])")
        max_diff_array = (refinement_limits[iref], refinement_limits[iref], refinement_limits[iref])
        SolidStateDetectors.refine!(sim, WeightingPotential, 1, max_diff_array)
        SolidStateDetectors.update_till_convergence!(sim, WeightingPotential, 1, depletion_handling=true, verbose=true)
        println("      r length: ",
            length(sim.weighting_potentials[1].grid.axes[1].ticks))
        println("      ϕ length: ",
            length(sim.weighting_potentials[1].grid.axes[2].ticks))
        println("      z length: ",
            length(sim.weighting_potentials[1].grid.axes[3].ticks))
    end
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


#  ============ PLOT ============

# Coordinate rettangolo
x_min, x_max = 0.015, 0.020
y_min, y_max = 0.02, 0.03
x_rect = [x_min, x_max, x_max, x_min, x_min]
y_rect = [y_min, y_min, y_max, y_max, y_min]

p = plot(sim.weighting_potentials[1],
    contours_equal_potential=true,
    linecolor=:white,
    levels=5,
    title="max tick = $max_tick_fin"
)

plot!(sim.detector,
    st=:slice,
    φ=90u"°",
    legend=false
)

plot!(x_rect, y_rect,
    seriestype=:shape,
    linecolor=:white,
    lw=1.5,
    fillalpha=0,
    label=""
)




savefig(p, "notebooks/plots/new_max_tick_function/$max_tick_fin.png")

#=

n = length(max_tick_fin_array)

ncols = 4
nrows = ceil(Int, n / ncols)
plt = plot(layout=(nrows, ncols), size=(1200, 800 * nrows))

# Coordinate rettangolo
x_min, x_max = 0.015, 0.020
y_min, y_max = 0.02, 0.04
x_rect = [x_min, x_max, x_max, x_min, x_min]
y_rect = [y_min, y_min, y_max, y_max, y_min]




for (i, max_tick_fin) in enumerate(max_tick_fin_array)
    println("⚡ Upoload simulation saved in : $save_sim_path")
    sim = ssd_read(save_sim_path, Simulation)
    println("===================================================")
    println("Initializing the state with a default grid...")
    new_grid = Grid(sim)
    apply_initial_state!(sim, WeightingPotential, 1, new_grid) # optional
    SolidStateDetectors.update_till_convergence!(sim, WeightingPotential, 1, depletion_handling=true, verbose=true)
    println("Initial grid initial state: ", length(sim.weighting_potentials[1].grid.axes[1].ticks))

    println("===================================================")
    println("max_tick_fin = $max_tick_fin")

    # Costruzione gerarchica dei max_tick
    max_tick_array = (
        max_tick_fin * fraction[1] * fraction[2] * fraction[3],
        max_tick_fin * fraction[1] * fraction[2],
        max_tick_fin * fraction[1],
        max_tick_fin
    )

    for max_tick in max_tick_array
        println("  Applying max tick distance = $max_tick")
        # qui faccio il refinemt sulla griglia
        SolidStateDetectors.refine_max_tick!(
            sim,
            WeightingPotential,
            1,
            (max_tick, max_tick, max_tick),
            (1e-15u"mm", 1e-15u"mm", 1e-15u"mm")
        )

        SolidStateDetectors.update_till_convergence!(
            sim,
            WeightingPotential,
            1,
            depletion_handling=true,
            verbose=true
        )

        n_refinement_steps = length(refinement_limits)
        # qui faccio il refinement sul potenziale (refinement classico)
        for iref in 1:n_refinement_steps
            max_diff_array = ntuple(_ -> refinement_limits[iref], 3)
            SolidStateDetectors.refine!(sim, WeightingPotential, 1, max_diff_array)

            SolidStateDetectors.update_till_convergence!(
                sim,
                WeightingPotential,
                1,
                depletion_handling=true,
                verbose=true
            )

            println("      r length: ",
                length(sim.weighting_potentials[1].grid.axes[1].ticks))
            println("      ϕ length: ",
                length(sim.weighting_potentials[1].grid.axes[2].ticks))
            println("      z length: ",
                length(sim.weighting_potentials[1].grid.axes[3].ticks))
        end
    end
    wp = sim.weighting_potentials[1]

    r_ticks = wp.grid.axes[1].ticks   # r
    z_ticks = wp.grid.axes[3].ticks   # z

    rmin, rmax = ustrip.(u"m", (0.015u"m", 0.02u"m"))
    zmin, zmax = ustrip.(u"m", (0.02u"m", 0.04u"m"))

    r_inds = findall(r -> rmin ≤ r ≤ rmax, r_ticks)
    z_inds = findall(z -> zmin ≤ z ≤ zmax, z_ticks)

    isempty(r_inds) && error("No points found in the r interval")
    isempty(z_inds) && error("No points found in the z interval")

    vals = wp.data[r_inds, 1, z_inds]
    min_wp = minimum(vals)

    println(">>> Min WeightingPotential in the test volume = $min_wp")

    # === UN SOLO PLOT PER max_tick_fin ===
    plot!(plt[i], sim.weighting_potentials[1],
        contours_equal_potential=true,
        linecolor=:white,
        levels=5,
        title="max tick fin = $max_tick_fin"
    )

    plot!(plt[i], sim.detector,
        st=:slice,
        φ=90u"°",
        legend=false
    )

    plot!(plt[i], x_rect, y_rect,
        seriestype=:shape,
        linecolor=:white,
        lw=1.5,
        fillalpha=0,
        label=""
    )
end

savefig(plt, "weighting_potential_vs_max_tick_fin.png")
=#