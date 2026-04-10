# =========================
# 📦 IMPORT PACKAGES
# =========================
using Pkg
Pkg.activate(".")

using Plots
using Unitful
using SolidStateDetectors
using JSON
using LegendHDF5IO
using ProgressMeter
using HDF5
using SpecialFunctions
using SolidStateDetectors: Electron, Hole

const T = Float64

println(pathof(SolidStateDetectors))

# =========================
# ⚙️ CONFIGURAZIONE
# =========================
det_z = det_r = 10u"mm"
z_draw = det_z / 2
pn_r = 8.957282u"mm" # this one was calculated by searching the zero impurity point (displayed in the following section)

recalculate = false
save_sim_path = "sim-saved/sim.h5"


# =========================
# 🔄 LOAD / CREATE SIMULATION
# =========================
if isfile(save_sim_path) && !recalculate
    println("⚡ Loading saved simulation...")
    sim = ssd_read(save_sim_path, Simulation)
else
    println("🔧 Creating new simulation...")

    sim = Simulation{T}(SSD_examples[:TrueCoaxial])
    cdm = sim.detector.semiconductor.charge_drift_model
    depth_list = 0u"mm":0.01u"mm":(det_r-pn_r)
    mobility_list = map(depth -> let pt::CartesianPoint{T} = CartesianPoint(det_r - depth, 0, z_draw)
            (µe=SolidStateDetectors.calculate_mobility(cdm, pt, Electron) * 10000u"cm^2/(V*s)",
                µh=SolidStateDetectors.calculate_mobility(cdm, pt, Hole) * 10000u"cm^2/(V*s)")
        end, depth_list)


    # calcolo il potenziale elettrico nella griglia rozza di default
    calculate_electric_potential!(sim, max_n_iterations=10, grid=Grid(sim), verbose=false, depletion_handling=true)

    g = sim.electric_potential.grid
    ax1, ax2, ax3 = g.axes

    bulk_tick_dis = 0.05u"mm"
    dl_tick_dis = 0.01u"mm"

    user_additional_ticks_ax1 = sort(vcat(ax1.interval.left*u"m":bulk_tick_dis:pn_r, pn_r:dl_tick_dis:ax1.interval.right*u"m"))
    user_ax1 = typeof(ax1)(ax1.interval, SolidStateDetectors.to_internal_units.(user_additional_ticks_ax1))
    user_g = typeof(g)((user_ax1, ax2, ax3))
    calculate_electric_potential!(sim, refinement_limits=0.1, grid=user_g, depletion_handling=true)
    calculate_electric_field!(sim)
    calculate_weighting_potential!(sim, 1, depletion_handling=true)
    calculate_weighting_potential!(sim, 2, depletion_handling=true)

    println("💾 Saving simulation...")
    ssd_write(save_sim_path, sim)
end

cdm = sim.detector.semiconductor.charge_drift_model

τh_fixed = sim.detector.semiconductor.charge_trapping_model.inactive_charge_trapping_model.τh .* 1e9u"ns"

#=
plot(
    begin
        imp = plot(sim.imp_scale, φ=0, xunit=u"mm", yunit=u"mm", title="impurity scale")
        vline!([pn_r], lw=2, ls=:dash, color=:darkred, label="PN junction boundary", legendfontsize=6)
    end,
    begin
        plot(sim.point_types, φ=0, xunit=u"mm", yunit=u"mm", title="point types")
        vline!([pn_r], lw=2, ls=:dash, color=:darkred, label="PN junction boundary", legendfontsize=6)
    end,
    begin
        plot(sim.electric_potential, xunit=u"mm", yunit=u"mm", title="electric potential")
        vline!([pn_r], lw=2, ls=:dash, color=:darkred, label="PN junction boundary", legendfontsize=6)
    end,
    begin
        plot(sim.electric_field, xunit=u"mm", yunit=u"mm", title="electric field", clims=(0, 100 * 2000))
        vline!([pn_r], lw=2, ls=:dash, color=:darkred, label="PN junction boundary", legendfontsize=6)
    end,
    size=(800, 600), layout=(2, 2),
)=#



# ===
# OLD CCE curves
# =======

# FONDAMENTALE 
dx = 0.1u"mm"
depth_list = 0.1u"mm":dx:(det_r-pn_r)
println("Simulating CCE with default τ")
totTime = 10000u"ns"
totEnergy = 5u"keV" # --> simulating ~339 carrier pairs
N = Int(totEnergy ÷ 2.95u"eV")  # 2.95 keV is the average energy needed to create an electron-hole pair in germanium


cdm = sim.detector.semiconductor.charge_drift_model
pbar = Progress(length(depth_list))
pulse_plot = plot()
eff_list_costτ = map(depth -> begin
        r = det_r - depth
        energy_depos = fill(2.95u"eV", N)
        starting_positions = repeat([CartesianPoint(r, 0, z_draw)], N)
        evt = Event(starting_positions, energy_depos)
        simulate!(evt, sim, Δt=1u"ns", max_nsteps=round(Int, totTime / 1u"ns"), diffusion=true, end_drift_when_no_field=false, self_repulsion=false)
        pulse = evt.waveforms[1]
        plot!(pulse_plot, pulse, label="Depth: $(round(typeof(depth), depth, digits = 1))", lw=2, yunit=u"eV/V")
        next!(pbar)
        maximum(pulse.signal) / N

    end, depth_list)

eff_list_tau_cost = ustrip.(eff_list_costτ)  # rimuove l'unità eV
eff_list_tau_cost = sqrt.(eff_list_tau_cost ./ N)

plot!(pulse_plot, legend=:topright, xlabel="Time / ns", ylabel="Amplitude / e", unitformat=:nounit)
cce_plot = plot(depth_list, eff_list_costτ, xlabel="Depth to surface / mm", ylabel="Charge collection efficiency", lw=2, color=:black, label="", unitformat=:nounit)
plot(pulse_plot, cce_plot, layout=(1, 2), size=(1000, 400), margin=5Plots.mm)


# ======
# Mobility curves
# =====
println("Saving mobility curves")
mobility_list = map(depth -> let pt::CartesianPoint{T} = CartesianPoint(det_r - depth, 0, z_draw)
        (µe=SolidStateDetectors.calculate_mobility(cdm, pt, Electron) * 10000u"cm^2/(V*s)",
            µh=SolidStateDetectors.calculate_mobility(cdm, pt, Hole) * 10000u"cm^2/(V*s)")
    end, depth_list);

µe_list = [m.µe for m in mobility_list]   # lista delle µe
µh_list = [m.µh for m in mobility_list]   # lista delle µh



# ===========
# Diffusion characteristics
# ===========


Δt = 1e-9u"s"

T_ann = sim.detector.semiconductor.impurity_density_model.surface_imp_model.lithium_annealing_temperature * u"K"  # K
t_ann = sim.detector.semiconductor.impurity_density_model.surface_imp_model.lithium_annealing_time * u"s"           # s

T_diff = sim.detector.semiconductor.temperature * u"K"  # temperatura del rivelatore in K
q = 1.602e-19u"C"    # carica elementare
kB = 1.380649e-23u"J/K"
#=
Dh = µh_list .* (kB * T_diff / q)
De = µe_list .* (kB * T_diff / q)
println("V * C / J = 1")

V_factor = (kB * T_diff / q) |> u"V"

# Calcolo della diffusività in cm^2/s
Dh = µh_list .* V_factor  #  cm^2/s
De = µe_list .* V_factor

σ_h = sqrt.(6 .* Dh .* Δt)
σ_e = sqrt.(6 .* De .* Δt)  # cm
=#

# =========
# Lithium concentration profile
# =========
println("Calculating the lifetime profile τ(x)")


R = 1.98u"cal/(mol*K)"
H = 11800u"cal/mol"
D0 = 2.5e-1u"mm^2/s"
Ns0 = 10^(21.27) * u"cm^-3"
Ns = Ns0 * 10^(-2610u"K" / T_ann)
D_Li = D0 * exp(-H / (R * T_ann))

r_Li = uconvert(u"m", 0.002u"cm")


# ==========
# carriers lifetime definition
# ==========
function τ_hole(x, α=1.6e-11)
    # --- Constants (SI) ---
    m0 = 9.11e-31u"kg"
    m_eff_hole = 0.21 * m0
    kB = 1.38e-23u"J/K"
    # Convert radius to meters (IMPORTANT)
    r_Li = uconvert(u"m", 0.002u"cm")
    # --- Thermal velocity (m/s) ---
    v_th = sqrt(3 * kB * T_diff / m_eff_hole)
    # --- Cross section (m²) ---
    σ_trap = π * r_Li^2
    # --- Dopant concentration (force to m⁻³) ---
    Nd_val = Ns .* erfc.(x ./ (2 * sqrt(D_Li * t_ann)))
    Nd_val = uconvert.(u"m^-3", Nd_val)
    # --- Lifetime (seconds) ---
    τ = 1.0 ./ (α .* Nd_val .* v_th .* σ_trap)
    return uconvert.(u"s", τ)
end
function τ_electron(x, α=1.6e-11)
    # --- Constants (SI) ---
    m0 = 9.11e-31u"kg"
    m_eff_electron = 0.12 * m0
    kB = 1.38e-23u"J/K"
    # Convert radius to meters (IMPORTANT)
    r_Li = uconvert(u"m", 0.002u"cm")
    # --- Thermal velocity (m/s) ---
    v_th = sqrt(3 * kB * T_diff / m_eff_electron)
    # --- Cross section (m²) ---
    σ_trap = π * r_Li^2
    # --- Dopant concentration (force to m⁻³) ---
    Nd_val = Ns .* erfc.(x ./ (2 * sqrt(D_Li * t_ann)))
    Nd_val = uconvert.(u"m^-3", Nd_val)
    # --- Lifetime (seconds) ---
    τ = 1.0 ./ (α .* Nd_val .* v_th .* σ_trap)
    return uconvert.(u"s", τ)
end


# =========
# New CCE curves with a position-dependent lifetime
# =========
cdm = sim.detector.semiconductor.charge_drift_model
println("Simulating CCE with τ(x)")
pulse_plot = plot()
τ_listh = Float64[]
τ_liste = Float64[]
pbar = Progress(length(depth_list))
eff_list_taux = map(depth -> begin
        # Calcolo della posizione di partenza
        r = det_r - depth

        # τh per questa profondità
        τh = ustrip(u"s", τ_hole(depth))
        τe = ustrip(u"s", τ_electron(depth))
        old_ctm = sim.detector.semiconductor.charge_trapping_model

        new_model = ConstantLifetimeChargeTrappingModel(
            τh,
            τe
        )

        new_ctm = CombinedChargeTrappingModel(
            old_ctm.bulk_charge_trapping_model,
            new_model,
            old_ctm.inactive_layer_geometry
        )

        sim.detector = SolidStateDetector(sim.detector, new_ctm)

        println("τh = ", sim.detector.semiconductor.charge_trapping_model.inactive_charge_trapping_model.τh)
        println("τe = ", sim.detector.semiconductor.charge_trapping_model.inactive_charge_trapping_model.τe)
        push!(τ_listh, sim.detector.semiconductor.charge_trapping_model.inactive_charge_trapping_model.τh)
        push!(τ_liste, sim.detector.semiconductor.charge_trapping_model.inactive_charge_trapping_model.τe)
        # Creazione eventi
        energy_depos = fill(2.95u"eV", N)
        starting_positions = repeat([CartesianPoint(r, 0, z_draw)], N)
        evt = Event(starting_positions, energy_depos)

        # Simulazione
        simulate!(evt, sim, Δt=1u"ns",
            max_nsteps=round(Int, totTime / 1u"ns"),
            diffusion=true,
            end_drift_when_no_field=false,
            self_repulsion=false)

        # Estrazione del pulse
        pulse = evt.waveforms[1]
        depth_val = ustrip(u"mm", depth)
        # Plot del pulse
        plot!(pulse_plot, pulse, label="Depth: $(round(depth_val, digits=1))", lw=2, yunit=u"eV/V")
        next!(pbar)
        # Calcolo efficienza di raccolta
        maximum(pulse.signal) / N

    end, depth_list)
# efficienza senza unità
eff_list_taux_num = ustrip.(eff_list_taux)  # rimuove l'unità eV
err_eff_taux = sqrt.(eff_list_taux_num ./ N)
#err_eff_taux = sqrt.(ustrip.(u"eV", eff_list_taux) ./ N)
print("miao")

plot!(pulse_plot, legend=:topright, xlabel="Time / ns", ylabel="Amplitude / e", unitformat=:nounit)
cce_plot_taux = plot(depth_list, eff_list_taux, xlabel="Depth to surface / mm", ylabel="Charge collection efficiency", lw=2, color=:black, label="", unitformat=:nounit)
plot(pulse_plot, cce_plot_taux, layout=(1, 2), size=(1000, 400), margin=5Plots.mm)



# confronto con risultato di Claudia

plot(
    depth_list,
    eff_list_costτ,
    yerr=eff_list_tau_cost,
    lw=2,
    label="τ = $τh_fixed",
    unitformat=:nounit,
    xticks=0:0.1:maximum(ustrip.(depth_list)),  # tick x ogni 1 mm
    yticks=0:0.1:1,                           # tick y ogni 0.05
    #elw=0.5,                         # barre di errore sottili
    xlabel="Depth to surface / mm",
    ylabel="Charge collection efficiency",
    legend=:topleft,
    frame=:box,)

plot!(
    depth_list,
    eff_list_taux,
    yerr=err_eff_taux,
    lw=2,
    label="τ(x)",
    unitformat=:nounit,
    #elw=0.5,
)