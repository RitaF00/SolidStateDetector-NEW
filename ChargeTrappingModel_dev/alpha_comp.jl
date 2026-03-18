using Random
using Plots
using SpecialFunctions
using ProgressMeter
using Distributions
using Statistics

# -----------------------------
# Parametri fisici e griglia
# -----------------------------
Lx = 0.2      # cm
Ly = 0.2
Lz = 0.2

dx = 0.001   # slice per diffusione Li (1 µm)

FCCD_cm = 0.1
PN_cm = 0.11

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
D = 28.9           # µm²/ns  --> va cambiato con la mobilità vera delle lacune
Δt = 1.0
t_max = 10000
Nt = Int(t_max / Δt)

σ = sqrt(2 * D * Δt) * 1e-4   # cm

# -----------------------------
# voxel precipitati
# -----------------------------
r_Li = 0.002  # 20 µm

# -----------------------------
# Funzioni
# -----------------------------
function generate_Li_grid(Lx, Ly, Lz, dx, α)
    nx = Int(Lx / r_Li)
    ny = Int(Ly / r_Li)
    nz = Int(Lz / r_Li)

    Li_grid = falses(nx, ny, nz)
    x_slices = 0:dx:Lx
    total_Li = 0

    for xi in x_slices
        Nd_val = Ns * erfc(xi / (2 * sqrt(D_Li * t_ann)))
        dV = dx * Ly * Lz
        λ = α * Nd_val * dV
        N_i = rand(Poisson(λ))
        total_Li += N_i

        for _ in 1:N_i
            x_pos = xi + rand() * dx
            y_pos = rand() * Ly
            z_pos = rand() * Lz
            ix = clamp(Int(floor(x_pos / r_Li)) + 1, 1, nx)
            iy = clamp(Int(floor(y_pos / r_Li)) + 1, 1, ny)
            iz = clamp(Int(floor(z_pos / r_Li)) + 1, 1, nz)
            Li_grid[ix, iy, iz] = true
        end
    end
    println("Numero totale di atomi Li generati = ", total_Li)
    return Li_grid, nx, ny, nz
end

function multiple_charges_trapping_3D(x_charges, N::Int=100, Li_grid=nothing)
    if isa(x_charges, Number)
        x_charges = fill(x_charges, N)
    end
    nx, ny, nz = size(Li_grid)
    collected_count = 0

    for n in 1:N
        x_charge = x_charges[n]
        y_charge = rand() * Ly
        z_charge = rand() * Lz

        for i in 1:Nt
            x_charge += σ * randn()
            y_charge += σ * randn()
            z_charge += σ * randn()

            x_charge = clamp(x_charge, 0, Lx)
            y_charge = clamp(y_charge, 0, Ly)
            z_charge = clamp(z_charge, 0, Lz)

            ix = clamp(Int(floor(x_charge / r_Li)) + 1, 1, nx)
            iy = clamp(Int(floor(y_charge / r_Li)) + 1, 1, ny)
            iz = clamp(Int(floor(z_charge / r_Li)) + 1, 1, nz)

            if Li_grid[ix, iy, iz]
                break
            end

            if x_charge >= FCCD_cm
                collected_count += 1
                break
            end
        end
    end
    return collected_count
end

# -----------------------------
# Variabili per alpha
# -----------------------------
alphas = [1e-8, 1e-9, 1e-10, 1e-11, 1e-12]
colors = [:deepskyblue, :slateblue2, :chocolate1, :chartreuse1, :orchid1]

x_pos = 0:0.001:0.15
N_charges = 1000
N_repeat = 7

# -----------------------------
# Loop su alpha e simulazione CCE
# -----------------------------
CCE_mean_list = []
CCE_std_list = []

for α in alphas
    println("Simulazione con α = ", α)
    Li_grid, nx, ny, nz = generate_Li_grid(Lx, Ly, Lz, dx, α)

    N_matrix = zeros(Int, length(x_pos), N_repeat)
    p = Progress(length(x_pos) * N_repeat)

    for (i, x0) in enumerate(x_pos)
        for j in 1:N_repeat
            N_matrix[i, j] = multiple_charges_trapping_3D(x0, N_charges, Li_grid)
            next!(p)
        end
    end

    CCE = N_matrix ./ N_charges
    push!(CCE_mean_list, vec(mean(CCE, dims=2)))
    push!(CCE_std_list, vec(std(CCE, dims=2)))
end

# -----------------------------
# Plot CCE
# -----------------------------
plt = plot(xlabel="depth (cm)", ylabel="CCE", frame=:box, size=(800, 600))

for i in 1:length(alphas)
    plot!(plt, x_pos, CCE_mean_list[i], lw=1.5, color=colors[i], label="α=$(alphas[i])")
    plot!(
        x_pos, CCE_mean_list[i],
        yerr=CCE_std_list[i],
        seriestype=:scatter,
        color=:black,
        marker=:circle,
        ms=1.5,
        label=false
    )
end

vline!(plt, [FCCD_cm], ls=:dash, color=:black, label="FCCD")
savefig(plt, "plot/CCE_vs_alpha.png")
display(plt)