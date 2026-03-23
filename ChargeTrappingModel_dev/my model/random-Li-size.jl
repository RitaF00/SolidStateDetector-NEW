using Random
using Plots
using SpecialFunctions
using ProgressMeter
using Distributions
using JSON
using StatsBase

# -----------------------------
# Parametri fisici
# -----------------------------
Lx = 0.2
Ly = 0.2
Lz = 0.2

dx = 0.0010   # slice per diffusione Li (1 µm)
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

α = 1.6e-11

# -----------------------------
# Diffusione carica
# -----------------------------
D = 28.9
Δt = 1.0
t_max = 10000
Nt = Int(t_max / Δt)

σ = sqrt(2 * D * Δt) * 1e-4  # cm

# -----------------------------
# Raggio precipitati
# -----------------------------
r_Li = 0.002  # 20 µm


function generate_Li_cells(Lx, Ly, Lz, dx, α)
    cell_size = r_Li
    nx = Int(Lx / cell_size)
    ny = Int(Ly / cell_size)
    nz = Int(Lz / cell_size)
    cells = [Vector{NTuple{4,Float64}}() for _ in 1:nx, _ in 1:ny, _ in 1:nz]

    x_slices = 0:dx:Lx
    total_Li = 0

    for xi in x_slices
        Nd_val = Ns * erfc(xi / (2 * sqrt(D_Li * t_ann)))
        dV = dx * Ly * Lz
        λ = α * Nd_val * dV
        N_i = rand(Poisson(λ))
        total_Li += N_i

        radii = [0.002, 0.006, 0.012]  # in cm
        weights = [0.7, 0.2, 0.1]

        for _ in 1:N_i
            x_pos = xi + rand() * dx
            y_pos = rand() * Ly
            z_pos = rand() * Lz
            r_Li = sample(radii, Weights(weights))

            ix = clamp(Int(floor(x_pos / cell_size)) + 1, 1, nx)
            iy = clamp(Int(floor(y_pos / cell_size)) + 1, 1, ny)
            iz = clamp(Int(floor(z_pos / cell_size)) + 1, 1, nz)

            push!(cells[ix, iy, iz], (x_pos, y_pos, z_pos, r_Li))
        end
    end

    println("α = $α: Numero totale di atomi Li generati = ", total_Li)
    return cells, nx, ny, nz, cell_size
end


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

            if trapped || x == 0
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

println("done")
cells, nx, ny, nz, cell_size = generate_Li_cells(Lx, Ly, Lz, dx, α)
n_precipitates = sum(length(cells[ix, iy, iz]) for ix in 1:nx, iy in 1:ny, iz in 1:nz)
println("Numero precipitati Li = ", n_precipitates)

xs_small, ys_small, zs_small = Float64[], Float64[], Float64[]
xs_med, ys_med, zs_med = Float64[], Float64[], Float64[]
xs_large, ys_large, zs_large = Float64[], Float64[], Float64[]

# cicla su tutte le celle
for ix in 1:nx, iy in 1:ny, iz in 1:nz
    for (x, y, z, r) in cells[ix, iy, iz]  # assumiamo tuple con raggio
        if r ≈ 0.002   # 20 μm
            push!(xs_small, x)
            push!(ys_small, y)
            push!(zs_small, z)
        elseif r ≈ 0.006  # 60 μm
            push!(xs_med, x)
            push!(ys_med, y)
            push!(zs_med, z)
        elseif r ≈ 0.012  # 120 μm
            push!(xs_large, x)
            push!(ys_large, y)
            push!(zs_large, z)
        end
    end
end

# scatter3d con colori diversi
p = scatter3d(xs_small, ys_small, zs_small,
    markersize=1.5,
    markercolor=:orange,
    label="r=20μm",
    camera=(25, 30),
    size=(800, 600)
)
scatter3d!(xs_med, ys_med, zs_med,
    markersize=2,
    markercolor=:magenta,
    label="r=60μm"
)
scatter3d!(xs_large, ys_large, zs_large,
    markersize=2.5,
    markercolor=:deepskyblue,
    label="r=120μm"
)

xlabel!("x (cm)")
ylabel!("y (cm)")
zlabel!("z (cm)")
xlims!(0, Lx)
ylims!(0, Ly)
zlims!(0, Lz)

xv = [0, Lx]
yv = [0, Ly]
zv = [0, Lz]

# Aggiungo le 12 linee degli spigoli
for x in xv, y in yv
    plot!(p, [x, x], [y, y], [0, Lz], linecolor=:black, label="")  # spigoli verticali
end
for x in xv, z in zv
    plot!(p, [x, x], [0, Ly], [z, z], linecolor=:black, label="")  # spigoli front/back
end
for y in yv, z in zv
    plot!(p, [0, Lx], [y, y], [z, z], linecolor=:black, label="")  # spigoli sinistra/destra
end

savefig(p, "plot/3D_Li_different_sizes.png")
display(p)



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




# ==== salvataggio del dizionario ====

filename = "random-Li.json"

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
