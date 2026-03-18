using Random
using Plots
using SpecialFunctions
using ProgressMeter
using Distributions
using JSON

# -----------------------------
# Parametri fisici e griglia
# -----------------------------
Lx = 0.2      # cm
Ly = 0.2
Lz = 0.2

dx = 0.001   # slice per diffusione Li (1 µm)
α = 2 * 1e-11     # thinning factor

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
# Generazione griglia Li
# -----------------------------
function generate_Li_grid(Lx, Ly, Lz, dx, α)

    nx = Int(Lx / r_Li)
    ny = Int(Ly / r_Li)
    nz = Int(Lz / r_Li)

    Li_grid = falses(nx, ny, nz)

    x_slices = 0:dx:Lx

    total_Li = 0   # contatore Li

    for xi in x_slices

        Nd_val = Ns * erfc(xi / (2 * sqrt(D_Li * t_ann)))

        dV = dx * Ly * Lz
        λ = α * Nd_val * dV

        N_i = rand(Poisson(λ))

        total_Li += N_i   # somma tutti i Li generati

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

    println("Numero totale di precipitati Li generati = ", total_Li)

    return Li_grid, nx, ny, nz
end


# -----------------------------
# Simulazione trapping cariche
# -----------------------------
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
            if x_charge == 0
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
# Generazione griglia Li
# -----------------------------
Li_grid, nx, ny, nz = generate_Li_grid(Lx, Ly, Lz, dx, α)
n_precipitates = count(Li_grid)

println("Numero precipitati Li = ", n_precipitates)
println("Matrice Li 3D generata!")

# -----------------------------
# Plot distribuzione Li
# -----------------------------
factor = 5

Li_plot = Li_grid

nxp, nyp, nzp = size(Li_plot)

xs = Float64[]
ys = Float64[]
zs = Float64[]

xs, ys, zs = Float64[], Float64[], Float64[]

for ix in 1:nx
    for iy in 1:ny
        for iz in 1:nz
            if Li_grid[ix, iy, iz]
                push!(xs, (ix - 0.5) * r_Li)
                push!(ys, (iy - 0.5) * r_Li)
                push!(zs, (iz - 0.5) * r_Li)
            end
        end
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
    zlim=(0, Lz), camera=(25, 30), size=(800, 800)
)
# Vertici del cubo
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
savefig(p, "3D_Li_rLi_$(round(r_Li*1e4))um.png")
display(p)
# -----------------------------
# Simulazione CCE
# -----------------------------
x_pos = 0:0.002:0.15

N_charges = 1000
N_repeat = 25

N_matrix = zeros(Int, length(x_pos), N_repeat)

p = Progress(length(x_pos) * N_repeat)

for (i, x0) in enumerate(x_pos)

    for j in 1:N_repeat

        N_matrix[i, j] = multiple_charges_trapping_3D(x0, N_charges, Li_grid)

        next!(p)

    end
end

CCE = N_matrix ./ N_charges
CCE_mean = vec(mean(CCE, dims=2))
CCE_std = vec(std(CCE, dims=2))

# -----------------------------
# Plot CCE
# -----------------------------
xticks_positions = 0:0.01:0.15
yticks_positions = 0:0.1:1
plt = plot(x_pos, CCE_mean, lw=5, color=:orange,
    xlabel="depth (cm)", ylabel="CCE", label="20 μm precipitates",
    xticks=(xticks_positions, string.(round.(xticks_positions, digits=2))),
    yticks=(yticks_positions, string.(round.(yticks_positions, digits=1))), frame=:box, size=(800, 600)
)
plot!(
    x_pos, CCE_mean,
    yerr=CCE_std,
    seriestype=:scatter,
    color=:black,
    marker=:circle,
    ms=1.5,
    label=false
)
vline!(plt, [FCCD_cm], ls=:dash, color=:green, label="FCCD")
savefig(plt, "plot/CCE_vs_precipitate_size.png")
display(plt)

#----dizionario---
results = Dict(
    "x_pos" => collect(x_pos),        # range → array
    "CCE_mean" => CCE_mean,
    "CCE_std" => CCE_std
)

# -----------------------------
# Salvataggio su file JSON
# -----------------------------
open("My-uncomplete-model.json", "w") do f
    JSON.print(f, results, 4)  # 4 = indentazione bella leggibile
end

println("File JSON salvato: My-uncomplete-model.json")