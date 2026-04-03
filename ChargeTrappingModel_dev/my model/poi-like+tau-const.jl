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
α = 1.6e-11   # valore di default

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
# Probabilità di trapping
# -----------------------------
function P_trapping(Δt)
    τ_ns = 800 #ns
    return 1 .- exp.(-Δt ./ τ_ns)
end

function generate_Li_cells(Lx, Ly, Lz, dx, α)

    cell_size = 0.0020 # creo una griglia di 20 μm in mm

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
        #println("Numero totale di atomi Li generati a $xi = ", N_i)

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

    println("Numero totale di atomi Li generati = ", total_Li)

    return cells, nx, ny, nz, cell_size
end

# -----------------------------
# Trapping cariche (con sfere)
# -----------------------------
function multiple_charges_trapping_3D(x_charges, N, cells, cell_size)

    if isa(x_charges, Number)
        x_charges = fill(x_charges, N)
    end

    nx, ny, nz = size(cells)
    collected_count = 0

    for n in 1:N

        x = x_charges[n]
        y = rand() * Ly
        z = rand() * Lz

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

            trapped = false

            for dxi in -1:1, dyi in -1:1, dzi in -1:1

                jx = ix + dxi
                jy = iy + dyi
                jz = iz + dzi

                if 1 ≤ jx ≤ nx && 1 ≤ jy ≤ ny && 1 ≤ jz ≤ nz

                    for (px, py, pz) in cells[jx, jy, jz]

                        if (x - px)^2 + (y - py)^2 + (z - pz)^2 < r_Li^2
                            if rand() < P_trapping(Δt)
                                trapped = true
                                break
                            end
                            break
                        end
                    end

                    if trapped
                        break
                    end
                end
            end

            if trapped
                break
            end

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


cells, nx, ny, nz, cell_size = generate_Li_cells(Lx, Ly, Lz, dx, α)
n_precipitates = sum(length(cells[ix, iy, iz]) for ix in 1:nx, iy in 1:ny, iz in 1:nz)
println("Numero precipitati Li = ", n_precipitates)

# -----------------------------
# Simulazione CCE
# -----------------------------
x_pos = 0:0.002:0.11

N_charges = 500
N_repeat = 25

N_matrix = zeros(Int, length(x_pos), N_repeat)

pbar = Progress(length(x_pos) * N_repeat)

for (i, x0) in enumerate(x_pos)
    for j in 1:N_repeat

        N_matrix[i, j] = multiple_charges_trapping_3D(x0, N_charges, cells, cell_size)

        next!(pbar)
    end
end

CCE = N_matrix ./ N_charges
CCE_mean = vec(mean(CCE, dims=2))
CCE_std = vec(std(CCE, dims=2))

# -----------------------------
# Plot CCE
# -----------------------------
plt = plot(x_pos .* 10, CCE_mean,
    lw=2,
    color=:orange,
    xlabel="depth (mm)",
    ylabel="CCE",
    label="20 μm precipitates",
    size=(450, 600),
    frame=:box
)

plot!(
    x_pos .* 10, CCE_mean,
    color=:orange,
    ribbon=CCE_std,
    fillalpha=0.32,
    ms=2,
    label=false
)

vline!(plt, [FCCD_cm .* 10], linestyle=:dash, label="FCCD")

#savefig(plt, "CCE.png")
display(plt)


filename = "json-CCE-file/Point-Li-and-const-tau-800ns.json"

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