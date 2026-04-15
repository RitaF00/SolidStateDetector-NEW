using Random
using Plots
using SpecialFunctions
using ProgressMeter
using Distributions
using JSON

# -----------------------------
# Parametri fisici
# -----------------------------
Lx = 0.2  #cm in mm
Ly = 0.2
Lz = 0.2

dx = 0.001   # slice per diffusione Li in cm (10 µm)
α = 1.6 * 1e-11     # thinning factor


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
D = 28.9 # μm^2/ns
Δt = 1.0  #ns
t_max = 10000

Nt = Int(t_max / Δt)
σ = sqrt(6 * D * Δt) * 1e-4  # cm

# -----------------------------
# Raggio precipitati
# -----------------------------
r_Li = 0.002  # 20 µm

# -----------------------------
# Generazione celle con Li puntiformi
# -----------------------------
"""
creo comunque una griglia per poter velocizzare il coedice e controlalree se la carica cade in un Li soltanto in celle adiacenti,
ma ogni  cella ha un nummro di Li diverso, e non fissato ad 1.
"""
function generate_Li_cells(Lx, Ly, Lz, dx, α)

    cell_size = 0.0020 # creo una griglia di 20 μm in cm

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
                            trapped = true
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

# -----------------------------
# Generazione Li
# -----------------------------
cells, nx, ny, nz, cell_size = generate_Li_cells(Lx, Ly, Lz, dx, α)
n_precipitates = sum(length(cells[ix, iy, iz]) for ix in 1:nx, iy in 1:ny, iz in 1:nz)
println("Numero precipitati Li = ", n_precipitates)

# -----------------------------
# Plot distribuzione Li
# -----------------------------
xs, ys, zs = Float64[], Float64[], Float64[]

for ix in 1:nx, iy in 1:ny, iz in 1:nz
    for (x, y, z) in cells[ix, iy, iz]
        push!(xs, x)
        push!(ys, y)
        push!(zs, z)
    end
end

p = scatter3d(xs, ys, zs,
    markersize=1,
    markercolor=:grey,
    xlabel="x (cm)",
    ylabel="y (cm)",
    zlabel="z (cm)",
    label="Li precipitates",
    xlim=(0, Lx),
    ylim=(0, Ly),
    zlim=(0, Lz),
    camera=(25, 30),
    size=(800, 800)
)

savefig(p, "3D_Li.png")
display(p)

# -----------------------------
# Simulazione CCE
# -----------------------------
x_pos = 0:0.001:0.11

N_charges = 500
N_repeat = 10

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
plt = plot(x_pos, CCE_mean,
    lw=2,
    color=:orange,
    xlabel="depth (cm)",
    ylabel="CCE",
    label="20 μm precipitates",
    size=(800, 600),
    frame=:box
)

plot!(
    x_pos, CCE_mean,
    yerr=CCE_std,
    seriestype=:scatter,
    ms=2,
    label=false
)

vline!(plt, [FCCD_cm], linestyle=:dash, label="FCCD")

savefig(plt, "CCE.png")
display(plt)



filename = "My-uncomplete-model-3D.json"

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
