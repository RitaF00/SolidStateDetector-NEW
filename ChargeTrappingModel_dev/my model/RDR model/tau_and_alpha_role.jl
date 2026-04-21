using Random
using Plots
using SpecialFunctions
using ProgressMeter
using Statistics
using Printf

gr()


"""
Nel codice:
- non utilizzo piu la distirbuzione di litio 
- modellizzo il materiale com eun mezzo efficce
- introduco il tempo di vita che dipende dal numero di trappole
- faccio variare α
- trovo le CCE
"""

# -----------------------------
# Parametri geometrici (cm)
# -----------------------------
Lx = 0.2
Ly = 0.2
Lz = 0.2
dx = 0.002
FCCD_cm = 0.1
x_RDR = 0.065

# -----------------------------
# Scan α
# -----------------------------
alpha_list = [1e-12, 1.6e-11, 1e-10, 1e-9]
colors = [:mediumpurple2, :deeppink, :deepskyblue, :orange]

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
# Probabilità trapping (con RDR)
# -----------------------------
function trapping_probability(x, Δt, α)

    if x >= x_RDR
        return 0.0
    end

    τ_ns = τ_hole(x, α) * 1e9

    if τ_ns <= 0
        return 1.0
    end

    return 1 - exp(-Δt / τ_ns)
end

# -----------------------------
# Simulazione diffusione
# -----------------------------
function multiple_charges_trapping_3D(x0, N, α)

    collected_count = 0

    for n in 1:N

        x = x0
        y = rand() * Ly
        z = rand() * Lz

        if x <= dx
            continue
        end

        for i in 1:Nt

            # diffusione
            x += σ * randn()
            y += σ * randn()
            z += σ * randn()

            x = clamp(x, 0, Lx)
            y = clamp(y, 0, Ly)
            z = clamp(z, 0, Lz)

            # trapping
            P = trapping_probability(x, Δt, α)

            if rand() < P
                break
            end

            # raccolta
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
x_pos = 0:dx:0.11
N_charges = 250
N_repeat = 25

CCE_all = Dict()

for (k, α) in enumerate(alpha_list)

    println("\nSimulazione α = ", α)

    CCE_matrix = zeros(length(x_pos), N_repeat)

    pbar = Progress(length(x_pos) * N_repeat)

    for (i, x0) in enumerate(x_pos)
        for j in 1:N_repeat

            collected = multiple_charges_trapping_3D(x0, N_charges, α)

            CCE_matrix[i, j] = collected
            next!(pbar)
        end
    end

    CCE = CCE_matrix ./ N_charges
    CCE_mean = vec(mean(CCE, dims=2))
    CCE_std = vec(std(CCE, dims=2))

    CCE_all[α] = (CCE_mean, CCE_std)
end

# -----------------------------
# PLOT
# -----------------------------
p = plot(
    xlabel="Depth (mm)",
    ylabel="CCE",
    frame=:box,
    size=(450, 600),
    dpi=300
)

for (k, α) in enumerate(alpha_list)

    if haskey(CCE_all, α)

        CCE_mean, CCE_std = CCE_all[α]

        plot!(
            p,
            x_pos .* 10,
            CCE_mean,
            lw=2,
            color=colors[k],
            ribbon=CCE_std,
            fillalpha=0.25,
            label=@sprintf("α = %.1e", α)
        )
    end
end

# linee di riferimento
vline!([FCCD_cm * 10], linestyle=:dot, color=:black, label="FCCD")
vline!([x_RDR * 10], linestyle=:dash, color=:red, label="RDR")

display(p)
savefig(p, "plot/CCE_vs_alpha_tau_model.png")