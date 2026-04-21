using Random
using SpecialFunctions
using Distributions
using Plots

# -----------------------------
# Parametri fisici
# -----------------------------
Lx = 0.2
Ly = 0.2
Lz = 0.2

α = 1.6e-11
r_Li = 0.002  # 20 µm

# Annealing Li
t_ann = 18 * 60
T_ann = 623
R = 1.98
H = 11800
D0 = 2.5e-3

Ns = 10^(21.27 - 2610 / T_ann)
D_Li = D0 * exp(-H / (R * T_ann))

# -----------------------------
# Funzione generazione precipitati
# -----------------------------
function generate_Li_cells(Lx, Ly, Lz, dx, α)
    cell_size = r_Li
    nx = Int(Lx / cell_size)
    ny = Int(Ly / cell_size)
    nz = Int(Lz / cell_size)
    cells = [Vector{NTuple{3,Float64}}() for _ in 1:nx, _ in 1:ny, _ in 1:nz]
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
            ix = clamp(Int(floor(x_pos / cell_size)) + 1, 1, nx)
            iy = clamp(Int(floor(y_pos / cell_size)) + 1, 1, ny)
            iz = clamp(Int(floor(z_pos / cell_size)) + 1, 1, nz)
            push!(cells[ix, iy, iz], (x_pos, y_pos, z_pos))
        end
    end
    return cells, total_Li
end

# -----------------------------
# Parametri dx
# -----------------------------
dx_values_cm = [0.0001, 0.0003, 0.0005, 0.0008, 0.001, 0.003, 0.005, 0.008, 0.01]  # cm dal più piccolo al più grande
dx_labels_um = ["1", "3", "5", "8", "10", "30", "50", "80", "100"]      # etichette µm corrispondenti
N_repeat = 50

# -----------------------------
# Simulazioni multiple per stimare sigma
# -----------------------------
n_prec_matrix = zeros(Int, length(dx_values_cm), N_repeat)

for (i, dx) in enumerate(dx_values_cm)
    for j in 1:N_repeat
        local cells, n_precipitates
        cells, n_precipitates = generate_Li_cells(Lx, Ly, Lz, dx, α)
        n_prec_matrix[i, j] = n_precipitates
    end
end

n_prec_mean = vec(mean(n_prec_matrix, dims=2))
n_prec_std = vec(std(n_prec_matrix, dims=2))

# -----------------------------
# Ordina dal più grande al più piccolo
# -----------------------------
sorted_indices = sortperm(dx_values_cm, rev=true)
dx_labels_sorted = dx_labels_um[sorted_indices]
n_prec_mean_sorted = n_prec_mean[sorted_indices]
n_prec_std_sorted = n_prec_std[sorted_indices]

# -----------------------------
# Plot equispaziato con etichette in µm
# -----------------------------
p = plot(1:length(dx_values_cm), n_prec_mean_sorted,
    #yerr=n_prec_std_sorted,
    #markerstrokewidth=1, linewidth=1,
    ribbon=n_prec_std_sorted,
    fillalpha=0.25,
    lw=2,
    marker=:circle,
    xlabel="dx (µm)",
    ylabel="Li precipitates",
    xticks=(1:length(dx_values_cm), dx_labels_sorted),
    legend=false,
    grid=true,
    frame=:box,
    size=(600, 400))

savefig(p, "plot/slice-comparison.png")
display(p)