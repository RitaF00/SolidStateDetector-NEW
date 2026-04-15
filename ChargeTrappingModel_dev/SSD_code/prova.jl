# =========================
# 📦 IMPORT
# =========================
using Pkg
Pkg.activate(".")

using Plots
using Unitful
using SolidStateDetectors

# =========================
# 📥 LOAD SIMULATION
# =========================
save_sim_path = "sim-saved/sim.h5"

sim = ssd_read(save_sim_path, Simulation)

# (opzionale ma utile se qualcosa manca)
# calculate_electric_potential!(sim)
# calculate_electric_field!(sim)
# calculate_weighting_potential!(sim)

# =========================
# 🎯 EVENT SETUP
# =========================
det_z = det_r = 10u"mm"
z_draw = det_z / 2
r = 9.5u"mm"
totEnergy = 1u"keV"
Epair = 2.95u"eV"

N = Int(totEnergy ÷ Epair)

starting_positions = repeat(
    [CartesianPoint(r, 0, z_draw)],
    N
)

energies = fill(Epair, N)

evt = Event(starting_positions, energies)

# =========================
# 🚀 DRIFT + DIFFUSION
# =========================
drift_charges!(evt, sim; diffusion=true)

# =========================
# 📊 PLOT TRAJECTORIES
# =========================
p = plot(title="Drift paths with diffusion", xlabel="x", ylabel="z")
scatter!([r], [z_draw], marker=:star, color=:black, label="initial position")
for dp in evt.drift_paths

    # elettroni
    plot!(p,
        [x.x .* 1000 for x in dp.e_path],
        [x.z .* 1000 for x in dp.e_path],
        label="",
        color=:deepskyblue,
        alpha=0.3
    )

    # buche
    plot!(p,
        [x.x .* 1000 for x in dp.h_path],
        [x.z .* 1000 for x in dp.h_path],
        label="",
        color=:orange,
        alpha=0.3
    )

end

display(p)