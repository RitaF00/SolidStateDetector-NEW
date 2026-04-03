using Random
using Plots
using SpecialFunctions
using ProgressMeter
using Distributions
using JSON

# -----------------------------
# Parametri fisici
# -----------------------------
Lx = 0.2  # mm
Ly = 0.2
Lz = 0.2

dx = 0.001   # slice per diffusione Li in mm (1 µm)
α = 1.6e-11     # thinning factor
FCCD_cm = 0.1
x_RDR = 0.065 # cm

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
D = 28.9       # diffusione
Δt = 1.0       # ns
t_max = 10000  # ns
Nt = Int(t_max / Δt)

σ = sqrt(6 * D * Δt) * 1e-4  # cm

# -----------------------------
# Raggio precipitati
# -----------------------------
r_Li = 0.002  # 20 μm in mm

# -----------------------------
# Funzione generazione celle Li
# -----------------------------
function generate_Li_cells(Lx, Ly, Lz, dx, α)
    # Griglia per celle di Li
    cell_size = 0.002 # 20 μm
    nx = Int(Lx / cell_size)
    ny = Int(Ly / cell_size)
    nz = Int(Lz / cell_size)

    cells = [Vector{NTuple{3,Float64}}() for _ in 1:nx, _ in 1:ny, _ in 1:nz]

    x_slices = 0:dx:Lx
    total_Li = 0

    for xi in x_slices
        # Genera Li solo se xi < x_RDR
        if xi < x_RDR
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
    end

    println("Numero totale di atomi Li generati = ", total_Li)
    return cells, nx, ny, nz, cell_size
end

# -----------------------------
# Funzione simulazione trapping 3D
# -----------------------------
function multiple_charges_trapping_3D(x_charges, N, cells, cell_size)
    if isa(x_charges, Number)
        x_charges = fill(x_charges, N)
    end

    nx, ny, nz = size(cells)
    collected_count = 0
    trapped_count = 0

    for n in 1:N
        x = x_charges[n]
        y = rand() * Ly
        z = rand() * Lz
        trapped = false

        for i in 1:Nt
            x += σ * randn()
            y += σ * randn()
            z += σ * randn()
            x = clamp(x, 0, Lx)
            y = clamp(y, 0, Ly)
            z = clamp(z, 0, Lz)

            ix = clamp(Int(floor(x / cell_size)) + 1, 1, nx)
            iy = clamp(Int(floor(y / cell_size)) + 1, 1, ny)
            iz = clamp(Int(floor(z / cell_size)) + 1, 1, nz)

            # controlla celle vicine
            for dxi in -1:1, dyi in -1:1, dzi in -1:1
                jx = ix + dxi
                jy = iy + dyi
                jz = iz + dzi
                if 1 ≤ jx ≤ nx && 1 ≤ jy ≤ ny && 1 ≤ jz ≤ nz
                    for (px, py, pz) in cells[jx, jy, jz]
                        if (x - px)^2 + (y - py)^2 + (z - pz)^2 < r_Li^2
                            trapped = true
                            break
                        end
                    end
                end
                if trapped
                    break
                end
            end

            if trapped
                trapped_count += 1
                break
            end

            if x >= FCCD_cm
                collected_count += 1
                break
            end
        end
    end

    return collected_count, trapped_count
end

# -----------------------------
# Generazione Li
# -----------------------------
cells, nx, ny, nz, cell_size = generate_Li_cells(Lx, Ly, Lz, dx, α)

# -----------------------------
# CCE e probabilità di trapping
# -----------------------------
dist = 0.002
x_pos = 0:dist:0.12
N_charges = 500
N_repeat = 25

CCE_mean = zeros(length(x_pos))
CCE_std = zeros(length(x_pos))
P_trap_mean = zeros(length(x_pos))

pbar = Progress(length(x_pos) * N_repeat)

for (i, x0) in enumerate(x_pos)
    CCE_vec = zeros(N_repeat)
    Ptrap_vec = zeros(N_repeat)

    for j in 1:N_repeat
        collected, trapped = multiple_charges_trapping_3D(x0, N_charges, cells, cell_size)
        CCE_vec[j] = collected / N_charges
        Ptrap_vec[j] = trapped / N_charges
        next!(pbar)
    end

    CCE_mean[i] = mean(CCE_vec)
    CCE_std[i] = std(CCE_vec)
    P_trap_mean[i] = mean(Ptrap_vec)
end

# -----------------------------
# Plot CCE + probabilità trapping
# -----------------------------
plt = plot(x_pos .* 10, CCE_mean,
    ribbon=CCE_std,
    lw=2,
    color=:orange,
    xlabel="depth (mm)",
    ylabel="CCE",
    label="CCE",
    frame=:box,
    size=(605, 400)
)

display(plt)

filename = "json-CCE-file/My-model_3D.json"
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