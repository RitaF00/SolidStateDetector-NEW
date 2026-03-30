using Random
using Plots
using SpecialFunctions
using ProgressMeter
using Distributions
using JSON
using Statistics
using LaTeXStrings

# -----------------------------
# Parametri geometrici (cm)
# -----------------------------
Lx = 0.2
Ly = 0.2
Lz = 0.2

dx = 0.001
α_global = 1.6e-11   # valore di default

FCCD_cm = 0.1

# -----------------------------
# Annealing Li
# -----------------------------
t_ann = 18 * 60
T_ann = 623

R = 1.98
H = 11800
D0 = 2.5e-3

Ns = 10^(21.27 - 2610 / T_ann)
D_Li = D0 * exp(-H / (R * T_ann))

# -----------------------------
# Diffusione carica
# -----------------------------
D = 28.9

Δt = 1.0  # ns (tempo step)
Nt = 10000

T_diff = 90 # K

σ = sqrt(2 * D * Δt) * 1e-4  # cm

# -----------------------------
# Precipitati (solo per σ_trap)
# -----------------------------
r_Li = 0.002  # cm

# -----------------------------
# Lifetime τ(x)
# -----------------------------
function τ_hole(x, α)

    m0 = 9.11e-31
    m_eff_hole = 0.3 * m0
    kB = 1.38e-23

    v_th = sqrt(3 * kB * T_diff / m_eff_hole) * 100

    σ_trap = π * r_Li^2

    Nd_val = Ns * erfc.(x ./ (2 * sqrt(D_Li * t_ann)))

    return 1.0 ./ (α .* Nd_val .* v_th .* σ_trap)
end

# -----------------------------
# Probabilità di trapping
# -----------------------------
function P_trapping(Δt, x, α)
    τ_ns = τ_hole(x, α) .* 1e9
    return 1 .- exp.(-Δt ./ τ_ns)
end

# -----------------------------
# Simulazione diffusione + trapping
# -----------------------------
function multiple_charges_trapping_3D(x_charges, N, α)

    if isa(x_charges, Number)
        x_charges = fill(x_charges, N)
    end

    collected_count = 0

    for n in 1:N

        x = x_charges[n]
        y = rand() * Ly
        z = rand() * Lz

        for i in 1:Nt

            # Diffusione
            x += σ * randn()
            y += σ * randn()
            z += σ * randn()

            x = clamp(x, 0, Lx)
            y = clamp(y, 0, Ly)
            z = clamp(z, 0, Lz)

            # Trapping continuo (τ-based)
            if rand() < P_trapping(Δt, x, α)
                break
            end

            # Boundary condizioni
            if x == 0
                break
            end

            if x >= FCCD_cm
                collected_count += 1
                break
            end
        end
    end

    return collected_count
end

# -----------------------------
# Simulazione CCE
# -----------------------------
alphas = [1e-8, 1e-9, 1e-10, 1e-11, 1e-12]
colors = [:deepskyblue, :royalblue1, :mediumpurple2, :deeppink, :orange]

alphas = [1.6 * 1e-11]
colors = [:orange]

x_pos = 0:0.002:0.11

N_charges = 500        # ↑ aumentato (più stabile)
N_repeat = 25          # ↑ aumentato

CCE_results = Dict()

plt = plot(
    xlabel="depth (mm)",
    ylabel="CCE",
    size=(450, 600),
    frame=:box,
    legend=:bottomright
)

for (k, α_val) in enumerate(alphas)

    println("\nSimulazione per α = $α_val")

    N_matrix = zeros(Int, length(x_pos), N_repeat)

    pbar = Progress(length(x_pos) * N_repeat)

    for (i, x0) in enumerate(x_pos)
        for j in 1:N_repeat

            N_matrix[i, j] = multiple_charges_trapping_3D(
                x0, N_charges, α_val
            )

            next!(pbar)
        end
    end

    CCE = N_matrix ./ N_charges
    CCE_mean = vec(mean(CCE, dims=2))
    CCE_std = vec(std(CCE, dims=2))

    # Salvataggio risultati
    CCE_results[string(α_val)] = Dict(
        "CCE_mean" => CCE_mean,
        "CCE_std" => CCE_std
    )

    # Linea media
    plot!(
        plt,
        x_pos .* 10,
        CCE_mean,
        lw=2,
        color=colors[k],
        label="α = $(α_val)"
    )

    # Banda errore
    plot!(
        plt,
        x_pos .* 10,
        CCE_mean,
        ribbon=CCE_std,
        fillalpha=0.25,
        color=colors[k],
        label=false
    )
end

# Linea FCCD
vline!(plt, [FCCD_cm .* 10], linestyle=:dash, label="FCCD")

display(plt)
savefig(plt, "plot/CCE_tau_model_apha=1,6e-11.png")

# -----------------------------
# (OPZIONALE) Salvataggio JSON
# -----------------------------
#=
filename = "tau_model_results.json"

results = Dict(
    "x_pos" => collect(x_pos),
    "CCE_results" => CCE_results
)

open(filename, "w") do f
    JSON.print(f, results, 4)
end

println("File salvato: $filename")
=#

# salvo solo l'ultimo valore di α
filename = "json-CCE-file/tau_model_results_1,6e-11.json"
results = Dict(
    "x_pos" => collect(x_pos),
    "FCCD_cm" => FCCD_cm,
    "alphas" => CCE_results
)

if isfile(filename)
    println("Il file $filename esiste già. Vuoi sovrascriverlo? (y/n)")
    risposta = readline()

    if lowercase(risposta) != "y"
        println("Operazione annullata, file non modificato.")
        return
    end
end

# Salvataggio
open(filename, "w") do f
    JSON.print(f, results, 4)
end

println("File JSON salvato: $filename")