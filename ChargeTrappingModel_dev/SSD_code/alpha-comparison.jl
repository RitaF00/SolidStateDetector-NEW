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

det_z = det_r = 10u"mm"
z_draw = det_z / 2
pn_r = 8.957282u"mm" # this one was calculated by searching the zero impurity point (displayed in the following section)

recalculate = false
save_sim_path = "sim-saved/sim.h5"

"""
In quetso codice vedo qual è l'effetto di alpha sulle curve CCE.
Per fare questo ho ripristinato la probailità di trapping non più come fattore Δt / τ
ma bbensì come 1 - exp ( Δt / τ)
"""

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


# ======================
# Definition of the simulation parameter and energy deposition
# ======================
dx = 0.05u"mm"
depth_list = 0.1u"mm":dx:(det_r-pn_r)
println("Simulating CCE with default τ")
totTime = 500u"ns"
totEnergy = 10u"keV" # --> simulating ~339 carrier pairs
N = Int(totEnergy ÷ 2.95u"eV")  # 2.95 keV is the average energy needed to create an electron-hole pair in germanium
cdm = sim.detector.semiconductor.charge_drift_model


# =========================
# Starting moility 
# =========================
println("Saving mobility curves")
mobility_list = map(depth -> let pt::CartesianPoint{T} = CartesianPoint(det_r - depth, 0, z_draw)
        (µe=SolidStateDetectors.calculate_mobility(cdm, pt, Electron) * 10000u"cm^2/(V*s)",
            µh=SolidStateDetectors.calculate_mobility(cdm, pt, Hole) * 10000u"cm^2/(V*s)")
    end, depth_list);

µe_list = [m.µe for m in mobility_list]   # lista delle µe
µh_list = [m.µh for m in mobility_list]   # lista delle µh




# ==============================
# Lithium concentration profile
# ==============================
println("Calculating the lifetime profile τ(x)")

Δt = 0.1e-9u"s"

T_ann = sim.detector.semiconductor.impurity_density_model.surface_imp_model.lithium_annealing_temperature * u"K"  # K
t_ann = sim.detector.semiconductor.impurity_density_model.surface_imp_model.lithium_annealing_time * u"s"           # s
T_diff = sim.detector.semiconductor.temperature * u"K"  # temperatura del rivelatore in K

R = 1.98u"cal/(mol*K)"
H = 11800u"cal/mol"
D0 = 2.5e-1u"mm^2/s"
Ns0 = 10^(21.27) * u"cm^-3"
Ns = Ns0 * 10^(-2610u"K" / T_ann)
D_Li = D0 * exp(-H / (R * T_ann))
# Lithium dimension
r_Li = uconvert(u"m", 0.002u"cm")

#=
println("Constant tau Li recombination sites")
kB = 1.38e-23u"J/K"
tau_cost = 1000 * 1e-9u"s"
N_sites = 1 / (tau_cost * sqrt(3 * kB * T_diff / (0.21 * 9.11e-31u"kg")) * π * r_Li^2) # m^-3
println(N_sites)
=#



# ==========
# carriers lifetime definition
# ==========
function τ_hole(x, α)
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
function τ_electron(x, α)
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
α_arr = [1e-6, 1e-7, 1e-8, 1e-10, 1e-12]

CCE_matrix = Vector{Vector{Float64}}()
CCE_err_matrix = Vector{Vector{Float64}}()

println("============= α curves ==============")


p = plot(depth_list, ustrip.(u"ns", τ_electron(depth_list, α_arr[1])), ylabel="τ(ns)", label="α = $(α_arr[1])", frame=:box, yscale=:log10)

for i in 2:length(α_arr)
    plot!(depth_list, ustrip.(u"ns", τ_electron(depth_list, α_arr[i])), label="α = $(α_arr[i])")
end
hline!([1000], label="constant lifetime", ls=:dash, color=:red)
display(p)


for α in α_arr
    println("Simulating CCE with α = ", α)
    #global cdm = sim.detector.semiconductor.charge_drift_model
    pbar = Progress(length(depth_list))
    eff_list = map(depth -> begin
            # Calcolo della posizione di partenza
            r = det_r - depth

            # τh per questa profondità
            τh = ustrip(u"s", τ_hole(depth, α))
            τe = ustrip(u"s", τ_electron(depth, α))
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
            next!(pbar)
            # Calcolo efficienza di raccolta
            maximum(pulse.signal) / N |> ustrip
        end, depth_list)
    push!(CCE_matrix, eff_list)
    err_eff = sqrt.(eff_list ./ N)
    push!(CCE_err_matrix, err_eff)
end

colors = [:deepskyblue, :royalblue1, :mediumpurple2, :deeppink, :orange]
#colors = [:deepskyblue, :mediumpurple2, :orange]

p = plot(
    xlabel="Depth to surface / mm",
    ylabel="Charge collection efficiency",
    frame=:box,
    legend=:topleft,
    title="Energy = $totEnergy",
    unitformat=:nounit
)

for (i, α) in enumerate(α_arr)
    plot!(
        p,
        depth_list,
        yerr=CCE_err_matrix[i],
        CCE_matrix[i],
        color=colors[i],
        lw=2,
        label="α = $(α)"
    )
end

display(p)
