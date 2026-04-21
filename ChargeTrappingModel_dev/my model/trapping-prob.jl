using Random
using Plots
using SpecialFunctions
using ProgressMeter
using Distributions
using JSON
using Statistics
using LaTeXStrings

"""
Codice in cui vedo come si comporta una probabilità di raccolta che dipende dallo spazio.
La carica (singola) viene generata e poi lasciata diffondere.
Ad ogni passo x, tau viene calcolata --> si calcola poi la prbabilità che la carica venga intrappolata o no.
Qui tau non è mai costante.
"""

# -----------------------------
# Parametri geometrici (cm)
# -----------------------------
Lx = 0.2
Ly = 0.2
Lz = 0.2
dx = 0.002
FCCD_cm = 0.1
x_RDR = 0.065  # 👉 regione senza trapping

# -----------------------------
# Parametro α
# -----------------------------
α_global = 1.6e-11

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
Δt = 1.0
Nt = 10000
T_diff = 90
σ = sqrt(2 * D * Δt) * 1e-4

# -----------------------------
# Precipitati Li
# -----------------------------
r_Li = 0.002

# -----------------------------
# Lifetime τ(x)
# -----------------------------
function τ_hole(x, α)
    m0 = 9.11e-31
    m_eff_hole = 0.21 * m0
    kB = 1.38e-23

    v_th = sqrt(3 * kB * T_diff / m_eff_hole) * 100
    σ_trap = π * r_Li^2
    Nd_val = Ns * erfc(x / (2 * sqrt(D_Li * t_ann)))

    return 1.0 / (α * Nd_val * v_th * σ_trap)
end

# -----------------------------
# Probabilità di trapping locale
# 👉 CON CONDIZIONE RDR
# -----------------------------
function trapping_probability(x, Δt, α)

    if x >= x_RDR
        return 0.0   # 👉 RDR: nessun trapping
    end

    τ_ns = τ_hole(x, α) * 1e9

    # sicurezza numerica
    if τ_ns <= 0
        return 1.0
    end

    return 1 - exp(-Δt / τ_ns)
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

        if x <= dx
            continue
        end

        for i in 1:Nt

            # Diffusione
            x += σ * randn()
            y += σ * randn()
            z += σ * randn()

            # Boundaries
            x = clamp(x, 0, Lx)
            y = clamp(y, 0, Ly)
            z = clamp(z, 0, Lz)

            # Probabilità di trapping
            P = trapping_probability(x, Δt, α)

            if rand() < P
                break
            end

            # Raccolta
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
α_val = α_global
println("\nSimulazione per α = $α_val")

x_pos = 0:dx:0.11
N_charges = 250
N_repeat = 30

N_matrix = zeros(Int, length(x_pos), N_repeat)
pbar = Progress(length(x_pos) * N_repeat)

for (i, x0) in enumerate(x_pos)
    for j in 1:N_repeat
        N_matrix[i, j] = multiple_charges_trapping_3D(x0, N_charges, α_val)
        next!(pbar)
    end
end

# -----------------------------
# Calcolo CCE
# -----------------------------
CCE = N_matrix ./ N_charges
CCE_mean = vec(mean(CCE, dims=2))
CCE_std = vec(std(CCE, dims=2))

# -----------------------------
# Probabilità di trapping modellata
# -----------------------------
P_trap_model = zeros(length(x_pos))

for (i, x) in enumerate(x_pos)
    if i == 1
        P_trap_model[i] = 1.0
    else
        P_trap_model[i] = trapping_probability(x, Δt, α_val)
    end
end

# -----------------------------
# Plot
# -----------------------------
bar_width = dx * 10

plt = plot(
    x_pos .* 10,
    CCE_mean,
    ribbon=CCE_std,
    fillalpha=0.25,
    color=:orange,
    label="CCE",
    gridalphan=0.3
)

plot!(
    x_pos .* 10,
    CCE_mean,
    lw=2,
    xlim=(0, 1.1),
    ylim=(0, 1.05),
    color=:orange,
    xlabel="Depth (mm)",
    ylabel="Trapping probability / CCE",
    label="",
    framestyle=:box,
    size=(450, 600)
)

bar!(
    x_pos .* 10,
    P_trap_model,
    lw=0,
    bar_width=bar_width,
    color=:grey,
    alpha=0.7,
    label="Trapping Probability"
)

# FCCD
vline!([FCCD_cm * 10], linestyle=:dot, color=:black, label="FCCD")

# 👉 RDR (utile visivamente)
vline!([x_RDR * 10], linestyle=:dash, color=:grey, label="RDR")

savefig(plt, "plot/CCE_and_trapping_probability_tau(x).png")
display(plt)

