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
dx = 0.002   # 20 µm in cm
FCCD_cm = 0.1
x_RDR = 0.065  # cm

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
Δt = 1.0   # ns
Nt = 10000
T_diff = 90  # K
σ = sqrt(2 * D * Δt) * 1e-4  # cm

# -----------------------------
# Precipitati Li
# -----------------------------
r_Li = 0.002  # cm

# -----------------------------
# Lifetime τ(x)
# -----------------------------
function τ_hole(x, α)
    m0 = 9.11e-31
    m_eff_hole = 0.21 * m0
    kB = 1.38e-23

    v_th = sqrt(3 * kB * T_diff / m_eff_hole) * 100  # cm/s
    σ_trap = π * r_Li^2
    Nd_val = Ns * erfc(x / (2 * sqrt(D_Li * t_ann)))

    return 1.0 / (α * Nd_val * v_th * σ_trap)
end

# -----------------------------
# Probabilità di trapping locale
# (modello puro, senza controllo su dx)
# -----------------------------
function trapping_probability(x, Δt, α)
    if x >= x_RDR
        return 0.0        # Regione attiva: nessun trapping
    else
        τ_ns = τ_hole(x, α) * 1e9
        return 1 - exp(-Δt / τ_ns)
    end
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

            # Raccolta carica
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
N_charges = 1000
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
# Probabilità di trapping modellata (solo primo bin = 1)
# -----------------------------
P_trap_model = zeros(length(x_pos))
for (i, x) in enumerate(x_pos)
    if i == 1
        P_trap_model[i] = 1.0        # solo superficie
    else
        P_trap_model[i] = trapping_probability(x, Δt, α_val)
    end
end

# -----------------------------
# Plot CCE + ribbon + istogramma trapping
# -----------------------------
bar_width = dx * 10

plt =# Ribbon errore CCE
    plot(
        x_pos .* 10,
        CCE_mean,
        ribbon=CCE_std,
        fillalpha=0.25,
        color=:orange,
        label=false
    )

# Ribbon errore CCE
plot!(
    x_pos .* 10,
    CCE_mean,
    ribbon=CCE_std,
    fillalpha=0.25,
    color=:orange,
    label=false
)

plot!(
    x_pos .* 10,
    CCE_mean,
    lw=2,
    xlim=(0, 1.2),
    ylim=(0, 1.05),
    color=:orange,
    xlabel="Depth (mm)",
    ylabel="Trapping probability / CCE",
    label="CCE",
    framestyle=:box,
    size=(450, 600)
)

# Istogramma probabilità di trapping (senza primo bin se vuoi)
bar!(
    x_pos .* 10,
    P_trap_model,
    lw=0,
    bar_width=bar_width,
    color=:grey,
    alpha=0.7,
    label="Trapping Probability"
)

# Linea FCCD
vline!([FCCD_cm * 10], linestyle=:dot, color=:black, label="FCCD")

savefig(plt, "plot/CCE_and_trapping_probability_tau(x).png")
display(plt)

# -----------------------------
# Salvataggio JSON
# -----------------------------
filename = "json-CCE-file/tau(x).json"
results = Dict(
    "x_pos" => collect(x_pos),        # range → array
    "CCE_mean" => CCE_mean,
    "CCE_std" => CCE_std
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

