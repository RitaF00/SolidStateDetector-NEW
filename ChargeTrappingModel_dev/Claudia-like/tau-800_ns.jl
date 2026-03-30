using Random
using Plots
using Statistics
using JSON
using ProgressMeter
using Distributions
# -----------------------------
# Parametri geometrici (cm)
# -----------------------------
Lx = 0.2
Ly = 0.2
Lz = 0.2

FCCD_cm = 0.1
T_diff = 90      # K
D = 28.9         # cm^2/s

Δt = 1.0         # ns
τ_const = 800.0  # ns
Nt = 100000       # number of diffusion steps
σ = sqrt(2 * D * Δt * 1e-9) * 1e-4  # cm, Δt in s




τ_mean = 800.0  # ns
τ_std = 80.0    # ns, 10% spread

tau_dist = Truncated(Normal(τ_mean, τ_std), 0, Inf)
# -----------------------------
# Probabilità di trapping costante
# -----------------------------


# -----------------------------
# Simulazione diffusione + trapping
# -----------------------------
function multiple_charges_trapping_3D_constant_tau(x_charges, N)
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

            # Clamp alle dimensioni
            x = clamp(x, 0, Lx)
            y = clamp(y, 0, Ly)
            z = clamp(z, 0, Lz)

            # Trapping costante
            τ_carica = rand(tau_dist)
            P_trap = 1 - exp(-Δt / τ_carica)
            if rand() < P_trap
                break
            end

            # Condizione di raccolta oltre FCCD
            if x >= FCCD_cm
                collected_count += 1
                break
            end

            # Boundary sinistro
            if x <= 0
                break
            end
        end
    end

    return collected_count
end

# -----------------------------
# Parametri simulazione
# -----------------------------
N_charges = 200
N_repeat = 5
x_pos = 0:0.002:0.11
pbar = Progress(length(x_pos) * N_repeat)


CCE_matrix = zeros(length(x_pos), N_repeat)

# -----------------------------
# Loop principale
# -----------------------------
for (i, x0) in enumerate(x_pos)
    for j in 1:N_repeat
        CCE_matrix[i, j] = multiple_charges_trapping_3D_constant_tau(x0, N_charges)
    end
end

CCE = CCE_matrix ./ N_charges
CCE_mean = vec(mean(CCE, dims=2))
CCE_std = vec(std(CCE, dims=2))

# -----------------------------
# Plot
# -----------------------------
x_mm = x_pos .* 10  # cm → mm

plt = plot(
    x_mm,
    CCE_mean,
    ribbon=CCE_std,
    xlabel="Depth [mm]",
    ylabel="CCE",
    lw=2,
    color=:deepskyblue,
    label="τ = 800 ns",
    fillalpha=0.3,
    frame=:box
)

vline!(plt, [FCCD_cm * 10], linestyle=:dash, label="FCCD")
display(plt)
savefig(plt, "plot/CCE_tau_constant.png")



filename = "constant-τ-800ns.json"

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