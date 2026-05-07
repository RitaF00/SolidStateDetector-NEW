
using Random
using Plots
using SpecialFunctions
using ProgressMeter
using Distributions
using JSON
using Unitful
using SolidStateDetectors
using SolidStateDetectors: Electron, Hole
T = Float64





# -----------------------------
# Parametri fisici e griglia
# -----------------------------
Lx = 0.2       # cm
Ly = 0.2
Lz = 0.2

r_Li = 0.002   # cm (20 μm)

FCCD_cm = 0.1
PN_cm = 0.11

# Parametri diffusione cariche
D = 2.89 * 1e-7          # μm^2/ns 2.89e-7 cm^2 ns^-1
Δt = 1.0           # ns
t_max = 10000      # ns
Nt = Int(t_max / Δt)
σ = sqrt(2 * D * Δt) * 1e-4  # cm (conversione μm -> cm)

max_concentration = 0.015   # 1.5 % per averte 1200 precipitati

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
# Funzione generazione Li 3D
# -----------------------------
function generate_Li_grid_precipitates(Lx_mm, Ly_mm, Lz_mm, precipitate_size_um)
    # -----------------------------
    # Parametri voxel e volume
    # -----------------------------
    voxel = precipitate_size_um         # lato del voxel = dimensione del precipitato (μm)
    nx = Int(Lx_mm * 1000 / voxel)     # numero voxel lungo x in μm
    ny = Int(Ly_mm * 1000 / voxel)     # numero voxel lungo y
    nz = Int(Lz_mm * 1000 / voxel)     # numero voxel lungo z

    # profondità del Li layer e concentrazione massima
    li_depth_um = 500.0                # 0.5 mm
    #max_concentration = 0.005           # 0.5 %


    # -----------------------------
    # Inizializzazione griglia
    # -----------------------------
    Li_grid = falses(nx, ny, nz)

    # -----------------------------
    # Popolamento griglia
    # -----------------------------
    for x in 1:nx
        depth = (x - 1) * voxel                                       # depth in μm --> cm
        conc = depth < li_depth_um ? max_concentration * erfc(depth / 10000 / (2 * sqrt(D_Li * t_ann))) : 0.0

        for y in 1:ny, z in 1:nz
            # 1 Li per voxel massimo
            if rand() < conc
                Li_grid[x, y, z] = true
            end
        end
    end

    return Li_grid, nx, ny, nz, voxel
end
# -----------------------------
# Funzione simulazione CCE
# -----------------------------
function multiple_charges_trapping_3D_linear(x_charges, N::Int=100, Li_grid=nothing, Lx=Lx, Ly=Ly, Lz=Lz)
    if isa(x_charges, Number)
        x_charges = fill(x_charges, N)
    elseif length(x_charges) != N
        error("lunghezza x_charges != N")
    end

    nx, ny, nz = size(Li_grid)
    r_Li = 0.002
    collected_count = 0

    for n in 1:N
        x_charge = x_charges[n]
        y_charge = rand() * Ly
        z_charge = rand() * Lz
        trapped = false

        if x_charge >= Lx
            collected_count += 1
            continue
        end

        for i in 1:Nt
            x_charge += σ * randn()
            y_charge += σ * randn()
            z_charge += σ * randn()

            x_charge = clamp(x_charge, 0, Lx) #clamp assicura che la carica rimanda nel volume
            y_charge = clamp(y_charge, 0, Ly)
            z_charge = clamp(z_charge, 0, Lz)

            ix = clamp(Int(floor(x_charge / r_Li)) + 1, 1, nx)  #la posizione diventa un indice della griglia
            iy = clamp(Int(floor(y_charge / r_Li)) + 1, 1, ny)
            iz = clamp(Int(floor(z_charge / r_Li)) + 1, 1, nz)

            if Li_grid[ix, iy, iz]
                trapped = true
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
# Generazione griglia Li 3D
# -----------------------------
Li_grid, nx, ny, nz, voxel = generate_Li_grid_precipitates(Lx * 10, Ly * 10, Lz * 10, r_Li * 1e4)
println("Griglia Li 3D generata!")

N_Li = count(Li_grid)
println("Numero totale di precipitati Li = ", N_Li)

# -----------------------------
# Estrazione coordinate per plot 3D
# -----------------------------
xs, ys, zs = Float64[], Float64[], Float64[]
for x in 1:nx, y in 1:ny, z in 1:nz
    if Li_grid[x, y, z]
        push!(xs, (x - 1) * voxel * 1e-4) # conversione μm -> cm
        push!(ys, (y - 1) * voxel * 1e-4)
        push!(zs, (z - 1) * voxel * 1e-4)
    end
end

# Plot 3D
p = scatter3d(xs, ys, zs, markersize=2, markercolor=:grey,
    xlabel="x (cm)", ylabel="y (cm)", zlabel="z (cm)",
    xlim=(0, Lx), ylim=(0, Ly), zlim=(0, Lz), grid=true,
    legend=false,
    camera=(25, 34), size=(800, 800),
    title="Number Li = $N_Li")

# Vertici del cubo
xv = [0, Lx]
yv = [0, Ly]
zv = [0, Lz]

# Aggiungo le 12 linee degli spigoli
for x in xv, y in yv
    plot!(p, [x, x], [y, y], [0, Lz], linecolor=:black)  # spigoli verticali
end
for x in xv, z in zv
    plot!(p, [x, x], [0, Ly], [z, z], linecolor=:black)  # spigoli front/back
end
for y in yv, z in zv
    plot!(p, [0, Lx], [y, y], [z, z], linecolor=:black)  # spigoli sinistra/destra
end
savefig(p, "plot/giov_$(max_concentration*100)%_erf.png")
display(p)


# -----------------------------
# Simulazione CCE
# -----------------------------
step = 0.002
x_pos = 0:step:0.11
N_charges = 500
N_repeat = 25


N_matrix = zeros(Int, length(x_pos), N_repeat)
p_bar = Progress(length(x_pos) * N_repeat, desc="Simulazione CCE")

for (i, x0) in enumerate(x_pos)
    for j in 1:N_repeat
        N_matrix[i, j] = multiple_charges_trapping_3D_linear(x0, N_charges, Li_grid)
        next!(p_bar)
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

vline!(plt, [0.05], ls=:dash, color=:black, label="Li layer")
savefig(plt, "plot/CCE_Giov_erf_$(max_concentration*100)%.png")
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
open("erf-$(max_concentration*100)%.json", "w") do f
    JSON.print(f, results, 4)  # 4 = indentazione bella leggibile
end

println("File JSON salvato: erf-$(max_concentration*100)%.json")
