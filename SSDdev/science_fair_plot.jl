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

T = Float32



save_sim_path = "saved_simulation/sim_sf.h5"

recalculate = true   # <<<<<<<< CAMBIA QUI
if isfile(save_sim_path) && !recalculate
    println("⚡ Upoload simulation saved in : $save_sim_path")
    sim_sf = ssd_read(save_sim_path, Simulation)

else

    println("🔧 New simulation for the electric potential...")
    sim_sf = Simulation{T}(SSD_examples[:InvertedCoax])
    sim_sf.detector = SolidStateDetector(sim_sf.detector, contact_id=2, contact_potential=3500u"V")
    grid = Grid(sim_sf, max_tick_distance=0.1u"mm")
    println("🔧 New simulation for the electric potential...")
    calculate_electric_potential!(sim_sf, grid=grid)
    println("🔧 New simulation for the electric field...")
    calculate_electric_field!(sim_sf)
    println("🔧 New simulation for the weighting potential...")
    #calculate_weighting_potential!(sim_sf, 1, refinement_limits=[0.2, 0.1]; grid=sim_sf.electric_potential.grid)
    #calculate_weighting_potential!(sim_sf, 2, refinement_limits=[0.2, 0.1]; grid=sim_sf.electric_potential.grid)

    println("💾 Saving simulation in $save_sim_path")
    ssd_write(save_sim_path, sim_sf)
end

sim = sim_sf






p = plot(
    ElectricPotential(
        sim.electric_potential.data .* replace(
            SolidStateDetectors.is_pn_junction_point_type.(
                length(sim.electric_potential.grid.φ) == 1 ?
                sim.point_types.data :
                SolidStateDetectors.get_2π_potential(sim.point_types).data
            ), 0 => T(NaN)
        ),
        sim.electric_potential.grid
    ),
    full_det=true,
    φ=0, color=:viridis,
    contours_equal_potential=true,
    linecolor=:White,
    legend=false,
    # tick verso l’interno
)


# Aggiungi slice del detector
plot!(sim.detector, full_det=true, st=:slice, lw=0.5, φ=0)



plot!(
    size=(2000, 1200),
    guidefontsize=12,
    title="",
    label="",
    tickfont=font(8),
    tick_length=6,         # lunghezza del tick in pixel
    tick_direction=:in)

x_ticks_cm = -0.04:0.02:0.04
y_ticks_cm = -0.010:0.02:0.09
x_labels = [@sprintf("%.2f", x) for x in x_ticks_cm]  # 1 decimale
y_labels = [@sprintf("%.2f", y) for y in y_ticks_cm]  # 2 decimali se vuoi


# Applica ticks distanziati
xticks!(p, x_ticks_cm, x_labels)
yticks!(p, y_ticks_cm, y_labels)
savefig(p, "science_fair/electric_pot.png")


p =
    plot(sim.electric_field, full_det=true, size=(500, 500), clims=(0, 5e5), legend=false)
plot_electric_fieldlines!(sim, full_det=true, sampling=1u"mm", offset=2u"mm")

# Aggiungi slice del detector
plot!(sim.detector, full_det=true, st=:slice, color=:viridis, lw=0.5, φ=0)
plot!(framestyle=:box,
    grid=true,
    size=(500, 500),
    guidefontsize=12,
    tickfont=font(8),
    tick_length=6,         # lunghezza del tick in pixel
    tick_direction=:in)

savefig(p, "science_fair/electric_field.png")


#=
p = plot(
    WeightingPotential(
        # multiply all entries in the electric potential that are NOT part of the semiconductor with NaN
        sim.weighting_potentials[1].data .* replace(SolidStateDetectors.is_pn_junction_point_type.(
                length(sim.weighting_potentials[1].grid.φ) == 1 ? sim.point_types.data : SolidStateDetectors.get_2π_potential(sim.point_types).data
            ), 0 => T(NaN)),
        sim.weighting_potentials[1].grid
    ), φ=0, contours_equal_potential=true, linecolor=:White, legend=false, framestyle=:box, grid=true
)
# Aggiungi slice del detector
plot!(sim.detector, st=:slice, lw=2, φ=0, xunit=u"m", yunit=u"m", zunit=u"m", xlabel="r ", ylabel="z ")
plot!(framestyle=:box,
    grid=true,
    gridalpha=5,
    size=(500, 500),
    guidefontsize=12,
    tickfont=font(8),
    tick_length=6,         # lunghezza del tick in pixel
    tick_direction=:in)

#x_ticks_cm = -0.04:0.02:0.04
#y_ticks_cm = -0.010:0.02:0.09
#x_labels = [@sprintf("%.2f", x) for x in x_ticks_cm]  # 1 decimale
#y_labels = [@sprintf("%.2f", y) for y in y_ticks_cm]  # 2 decimali se vuoi


# Applica ticks distanziati
#xticks!(p, x_ticks_cm, x_labels)
#yticks!(p, y_ticks_cm, y_labels)
savefig(p, "CM_napoli/weight_pot1.png")

p = plot(
    WeightingPotential(
        # multiply all entries in the electric potential that are NOT part of the semiconductor with NaN
        sim.weighting_potentials[2].data .* replace(SolidStateDetectors.is_pn_junction_point_type.(
                length(sim.weighting_potentials[2].grid.φ) == 1 ? sim.point_types.data : SolidStateDetectors.get_2π_potential(sim.point_types).data
            ), 0 => T(NaN)),
        sim.weighting_potentials[2].grid
    ), φ=0, contours_equal_potential=true, linecolor=:White, legend=false, framestyle=:box, grid=true
)
# Aggiungi slice del detector
plot!(sim.detector, st=:slice, lw=2, φ=0, xunit=u"m", yunit=u"m", zunit=u"m", xlabel="r ", ylabel="z ")
plot!(framestyle=:box,
    grid=true,
    gridalpha=5,
    size=(500, 500),
    guidefontsize=12,
    tickfont=font(8),
    tick_length=6,         # lunghezza del tick in pixel
    tick_direction=:in)

#x_ticks_cm = -0.04:0.02:0.04
#y_ticks_cm = -0.010:0.02:0.09
#x_labels = [@sprintf("%.2f", x) for x in x_ticks_cm]  # 1 decimale
#y_labels = [@sprintf("%.2f", y) for y in y_ticks_cm]  # 2 decimali se vuoi


# Applica ticks distanziati
#xticks!(p, x_ticks_cm, x_labels)
#yticks!(p, y_ticks_cm, y_labels)
savefig(p, "CM_napoli/weight_pot2.png")


starting_positions = [CartesianPoint{T}(0.020, 0, 0.015),
    CartesianPoint{T}(-0.015, -0.02, 0.045),
    CartesianPoint{T}(-0.022, 0.015, 0.025)]
energy_depos = T[1460, 609, 1000] * u"keV" # are needed later in the signal generation

evt = Event(starting_positions, energy_depos);
drift_charges!(evt, sim, Δt=5u"ns")
evt.drift_paths
p = plot(sim.detector, size=(600, 600),
    xunit=u"mm", yunit=u"mm", zunit=u"mm", xlims=(-40, 40), ylims=(-40, 40), zlims=(-10, 90), camera=(45, 15), legend=false
)
plot!(evt.drift_paths, linewidth=2, linestyle=:dash, legendfontsize=10)
plot!(vcat(evt.locations...), color=:orange, markersize=6, label="")
savefig(p, "science_fair/charge_drift.png")



simulate!(evt, sim, Δt=5u"ns")
p = plot(u"ns", u"fC", framestyle=:box,)
plot!(evt.waveforms, linewidth=4, legend=false)
savefig(p, "science_fair/charge_dcollection.png")
=#