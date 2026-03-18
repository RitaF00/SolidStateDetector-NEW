using Random
using Plots
using ProgressMeter
using Statistics

# -----------------------------
# Parametri fisici e griglia
# -----------------------------
Lx = 0.2       # cm, lunghezza totale
Ly = 0.2       # cm
Lz = 0.2       # cm

FCCD_cm = 0.1  # cm
li_depth_cm = 0.05  # 0.5 mm profondità layer Li

# Diffusione cariche
D = 28.9       # μm^2/ns
Δt = 1.0       # ns
t_max = 10000  # ns
Nt = Int(t_max / Δt)
σ = sqrt(2 * D * Δt) * 1e-4  # cm

# -----------------------------
# Funzione generazione Li 3D
# -----------------------------
function generate_Li_grid_precipitates(Lx, Ly, Lz, r_Li; C_tot=0.01, li_depth_cm=0.05)
    # numero voxel in ciascuna direzione
    nx = Int(ceil(Lx / r_Li))
    ny = Int(ceil(Ly / r_Li))
    nz = Int(ceil(Lz / r_Li))

    # voxel layer Li lungo x
    nx_li = Int(floor(li_depth_cm / r_Li))
    n_voxel_total = nx_li * ny * nz

    # numero voxel occupati da Li
    n_Li_voxel = round(Int, C_tot * n_voxel_total)

    # selezione casuale voxel
    indices = shuffle(1:n_voxel_total)
    Li_indices = indices[1:n_Li_voxel]

    # griglia vuota
    Li_grid = falses(nx, ny, nz)

    for idx in Li_indices
        ix = div(idx - 1, ny * nz) + 1
        rem1 = (idx - 1) % (ny * nz)
        iy = div(rem1, nz) + 1
        iz = rem1 % nz + 1
        Li_grid[ix, iy, iz] = true
    end
    return Li_grid, nx, ny, nz, r_Li
end

# -----------------------------
# Funzione simulazione CCE
# -----------------------------
function multiple_charges_trapping_3D_linear(x_charges, r_Li, N::Int=10, Li_grid=nothing, Lx=Lx, Ly=Ly, Lz=Lz)
    if isa(x_charges, Number)
        x_charges = fill(x_charges, N)
    elseif length(x_charges) != N
        error("lunghezza x_charges != N")
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
                break  # catturata
            end

            if x_charge >= FCCD_cm
                collected_count += 1
                break
            end
            if x_charge == 0
                break
            end
        end
    end


    return collected_count
end

# -----------------------------
# Parametri simulazione
# -----------------------------
step = 0.002          # 20 μm
x_pos = 0:step:0.1    # fino a 1 mm
N_charges = 1000       # 10 eventi per punto
N_repeat = 25          # ripetizioni


# Raggi Li: 20 μm, 60 μm, 120 μm
r_Li_values = [0.002, 0.006, 0.012]  # cm

CCE_results = Dict()
CCE_std_results = Dict()

# -----------------------------
# Loop principale
# -----------------------------
for r_Li in r_Li_values
    println("Simulazione per r_Li = $(r_Li*1e4) μm")

    # Generazione griglia Li per plot 3D
    Li_grid, nx, ny, nz, voxel = generate_Li_grid_precipitates(Lx, Ly, Lz, r_Li; C_tot=0.01, li_depth_cm=li_depth_cm)

    # Estrazione coordinate centro voxel
    xs, ys, zs = Float64[], Float64[], Float64[]
    for x in 1:nx, y in 1:ny, z in 1:nz
        if Li_grid[x, y, z]
            push!(xs, (x - 0.5) * voxel)
            push!(ys, (y - 0.5) * voxel)
            push!(zs, (z - 0.5) * voxel)
        end
    end

    # Plot 3D scatter
    p = scatter3d(xs, ys, zs, markersize=3, markercolor=:grey,
        xlabel="x (cm)", ylabel="y (cm)", zlabel="z (cm)",
        xlim=(0, Lx), ylim=(0, Ly), zlim=(0, Lz), legend=false,
        camera=(25, 34), size=(600, 600))
    savefig(p, "3D_Li_rLi_$(round(r_Li*1e4))um.png")
    display(p)

    # -----------------------------
    # Simulazione CCE
    # -----------------------------
    N_matrix = zeros(Int, length(x_pos), N_repeat)
    p_bar = Progress(length(x_pos) * N_repeat, desc="CCE r_Li=$(r_Li*1e4)μm")

    for (i, x0) in enumerate(x_pos)
        for j in 1:N_repeat
            # rigenera griglia casuale per ogni evento (come nell’articolo)
            Li_grid_event, _, _, _, _ = generate_Li_grid_precipitates(Lx, Ly, Lz, r_Li; C_tot=0.01, li_depth_cm=li_depth_cm)
            N_matrix[i, j] = multiple_charges_trapping_3D_linear(x0, r_Li, N_charges, Li_grid_event)
            next!(p_bar)
        end
    end

    CCE = N_matrix ./ N_charges
    CCE_results[r_Li] = vec(mean(CCE, dims=2))
    CCE_std_results[r_Li] = vec(std(CCE, dims=2))
end

# -----------------------------
# Plot finale CCE
# -----------------------------
colors = [:deepskyblue, :lawngreen, :deeppink]
xticks_positions = 0:0.01:0.15
yticks_positions = 0:0.1:1
plt = plot(
    xlabel="depth (cm)",
    ylabel="CCE",
    xticks=(xticks_positions, string.(round.(xticks_positions, digits=2))),
    yticks=(yticks_positions, string.(round.(yticks_positions, digits=1))),
    frame=:box,
    size=(800, 600)
)

for (i, r_Li) in enumerate(r_Li_values)
    plot!(x_pos, CCE_results[r_Li], lw=2.5, color=colors[i], label="$(round(r_Li*1e4)) μm")
    # barre di errore
    plot!(x_pos, CCE_results[r_Li], yerr=CCE_std_results[r_Li], seriestype=:scatter, marker=:circle, ms=2, color=:black, label=false)
end

vline!(plt, [li_depth_cm], ls=:dash, color=:black, label="Li layer depth")
savefig(plt, "CCE_vs_precipitate_size.png")
display(plt)