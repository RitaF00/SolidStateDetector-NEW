using Plots
using SpecialFunctions
using Distributions

gr()

# ============================================================
# SIMULAZIONE MONTE CARLO DELLA PRECIPITAZIONE DI PARTICELLE Li
# ============================================================
#
# Questo codice simula la formazione e distribuzione spaziale
# di particelle sferiche di litio all'interno di un volume 3D
# soggetto a un profilo di diffusione dipendente dalla profondità.
#
# Il modello combina:
# - Diffusione del litio (legge tipo erfc)
# - Nucleazione stocastica (distribuzione di Poisson)
# - Vincolo geometrico di non-sovrapposizione (hard-sphere packing)
#
# PROCEDURA PRINCIPALE:
# 1. Calcolo della profondità di saturazione del litio
#    in base al profilo diffusivo e a una soglia critica.
#
# 2. Generazione di particelle lungo la profondità:
#    - La densità locale dipende dalla concentrazione di Li
#    - Il numero di particelle segue una statistica di Poisson
#    - Il raggio delle particelle varia con la profondità
#
# 3. Inserimento delle particelle in un dominio 3D con:
#    - Subdivisione spaziale in celle (cell list)
#    - Controllo locale delle sovrapposizioni
#
# 4. Accettazione/rifiuto delle particelle in base a:
#    - Vincolo geometrico (no overlap)
#    - Limite massimo di packing volumetrico
#
# 5. Stima della frazione volumica tramite Monte Carlo:
#    - Campionamento casuale nello spazio
#    - Verifica appartenenza a particelle sferiche
#
# OUTPUT:
# - Distribuzione delle particelle accettate vs generate
# - Evoluzione con la profondità
# - Frazione volumica locale (%)
# - Confronto per diversi raggi delle particelle
#
# ============================================================

# ============================================================
# STRUTTURA PARTICELLA
# ============================================================
struct LiParticle
    x::Float64
    y::Float64
    z::Float64
    r::Float64
end

# ============================================================
# PARAMETRI
# ============================================================
Lx = 0.2
Ly = 0.2
Lz = 0.2

dx = 0.001
r_Li = 0.002

n_arr = [1, 1.5, 2, 3, 4, 6]

t_ann = 18 * 60
T_ann = 623

R = 1.98
H = 11800
D0 = 2.5e-3
Ns = 10^(21.27 - 2610 / T_ann)
D_Li = D0 * exp(-H / (R * T_ann))

Nd_saturation = 1e14
α = 1.6e-11





# ============================================================
# SATURATION DEPTH
# ============================================================
function compute_saturation_depth()
    depth = 0:dx:0.11
    Nd = Ns .* erfc.(depth ./ (2 * sqrt(D_Li * t_ann)))
    idx = findfirst(i -> Nd[i] <= Nd_saturation, eachindex(depth))
    return depth[idx]
end

x_sat = compute_saturation_depth()
xs_ref = collect(0:dx:x_sat)

# ============================================================
# RAGGIO
# ============================================================
r_Li_profile(x, n, sat) = x <= sat / 2 ? n * r_Li : r_Li

# ============================================================
# OVERLAP CHECK
# ============================================================
function is_overlapping(x, y, z, r, cells, nx, ny, nz, cell)

    ix = clamp(Int(floor(x / cell)) + 1, 1, nx)
    iy = clamp(Int(floor(y / cell)) + 1, 1, ny)
    iz = clamp(Int(floor(z / cell)) + 1, 1, nz)

    for dx in -1:1, dy in -1:1, dz in -1:1
        jx, jy, jz = ix + dx, iy + dy, iz + dz

        if 1 ≤ jx ≤ nx && 1 ≤ jy ≤ ny && 1 ≤ jz ≤ nz
            for p in cells[jx, jy, jz]
                if (x - p.x)^2 + (y - p.y)^2 + (z - p.z)^2 < (r + p.r)^2
                    return true
                end
            end
        end
    end

    return false
end

# ============================================================
# GENERAZIONE LITIO
# ============================================================
function generate_Li(n, sat)

    cell = 2 * n * r_Li

    nx = floor(Int, Lx / cell)
    ny = floor(Int, Ly / cell)
    nz = floor(Int, Lz / cell)

    cells = [Vector{LiParticle}() for _ in 1:nx, _ in 1:ny, _ in 1:nz]

    xs = collect(0:dx:x_sat)

    Ni_geom = zeros(Int, length(xs))
    Ni_accept = zeros(Int, length(xs))

    packing = 0.64

    for (i, xi) in enumerate(xs)

        if xi >= sat
            continue
        end

        Nd = Ns * erfc(xi / (2 * sqrt(D_Li * t_ann)))
        density = α * max(Nd - Nd_saturation, 0)

        λ = density * dx * Ly * Lz

        r = r_Li_profile(xi, n, sat)

        Vp = (4 / 3) * π * r^3
        Vs = dx * Ly * Lz

        Nmax = Int(floor(packing * Vs / Vp))
        Ntarget = min(rand(Poisson(λ)), Nmax)

        Ni_geom[i] = Ntarget

        accepted = 0
        trials = max(20 * Ntarget, 50)

        for _ in 1:trials
            accepted >= Ntarget && break

            x = xi + rand() * dx
            y = rand() * Ly
            z = rand() * Lz

            r = r_Li_profile(x, n, sat)

            if !is_overlapping(x, y, z, r, cells, nx, ny, nz, cell)

                ix = clamp(Int(floor(x / cell)) + 1, 1, nx)
                iy = clamp(Int(floor(y / cell)) + 1, 1, ny)
                iz = clamp(Int(floor(z / cell)) + 1, 1, nz)

                push!(cells[ix, iy, iz], LiParticle(x, y, z, r))
                accepted += 1
            end
        end

        Ni_accept[i] = accepted
    end

    return cells, Ni_geom, Ni_accept, xs, cell
end

# ============================================================
# MONTE CARLO FRAZIONE VOLUMICA
# ============================================================
function volume_fraction_MC(cells, xs, cell; Nsamp=3000)

    nx, ny, nz = size(cells)
    φ = zeros(length(xs))

    for (i, xi) in enumerate(xs)

        inside = 0

        for _ in 1:Nsamp

            x = xi + rand() * dx
            y = rand() * Ly
            z = rand() * Lz

            ix = clamp(Int(floor(x / cell)) + 1, 1, nx)
            iy = clamp(Int(floor(y / cell)) + 1, 1, ny)
            iz = clamp(Int(floor(z / cell)) + 1, 1, nz)

            found = false

            for dxi in -1:1, dyi in -1:1, dzi in -1:1
                jx, jy, jz = ix + dxi, iy + dyi, iz + dzi

                if 1 ≤ jx ≤ nx && 1 ≤ jy ≤ ny && 1 ≤ jz ≤ nz
                    for p in cells[jx, jy, jz]
                        if (x - p.x)^2 + (y - p.y)^2 + (z - p.z)^2 < p.r^2
                            inside += 1
                            found = true
                            break
                        end
                    end
                end

                found && break
            end
        end

        φ[i] = inside / Nsamp
    end

    return φ
end

# ============================================================
# RUN SIMULATION
# ============================================================
sat = compute_saturation_depth()
println("Saturation depth = ", sat * 10, " mm")

all_geom = []
all_acc = []
all_phi = []



for n in n_arr

    global xs_ref

    println(" ===============  Running for n = ", n, " =============== ")

    cells, Ng, Na, xs, cell = generate_Li(n, sat)



    if xs_ref === nothing
        xs_ref = xs
    end

    φ = volume_fraction_MC(cells, xs, cell)

    push!(all_geom, Ng)
    push!(all_acc, Na)
    push!(all_phi, φ)
end

x_mm = xs_ref .* 10
x_mm = xs_ref .* 10
x_split = sat * 10 / 2

# ============================================================
# PLOT 1: GRID 3x4
# ============================================================
l = @layout [grid(2, 3)]

p_grid = plot(layout=l, size=(1800, 750),
    margin=5Plots.mm,
    left_margin=10Plots.mm,
    bottom_margin=5Plots.mm,
    top_margin=5Plots.mm,
    frame=:box)
xticks_positions = 0:0.05:sat*10
for (i, n) in enumerate(n_arr)
    local yticks_positions
    yticks_positions = range(0, floor(Int, maximum(all_geom[i])), length=5)
    plot!(p_grid[i], x_mm, all_geom[i], color=:black, label="Requested (Pois + geom")
    plot!(p_grid[i], x_mm, all_acc[i],
        color=:orange,
        lw=2,
        xticks=xticks_positions,
        yticks=yticks_positions,
        gridalpha=0.05,
        ls=:dash, label="Accepted (no overlap)")
    title!(p_grid[i], "r=$(n * r_Li * 10000) μm")
    vline!(p_grid[i], [x_split], ls=:dash, color=:black, label="RDR/2")
    xlabel!(p_grid[i], "Depth (mm)")
    ylabel!(p_grid[i], "Number of particles")
end

display(p_grid)

# ============================================================
# PLOT 2: CONFRONTO GLOBALE
# ============================================================
l = @layout [grid(1, 2)]
p2 = plot(layout=l, size=(1000, 400),
    margin=5Plots.mm,
    left_margin=10Plots.mm,
    bottom_margin=5Plots.mm, frame=:box)
yticks_positions = 0:floor(Int, maximum(vcat(all_geom...)) / 15):maximum(vcat(all_geom...))
xticks_positions = 0:0.1:sat*10
for (i, n) in enumerate(n_arr)

    plot!(p2[1], x_mm, all_geom[i], xticks=xticks_positions, yticks=yticks_positions, gridalpha=0.05, label="r=$(n * r_Li * 10000) μm")
    plot!(p2[2], x_mm, all_acc[i], xticks=xticks_positions, yticks=yticks_positions, gridalpha=0.05, label="r=$(n * r_Li * 1000) μm")
    plot!(xlabel="Depth (mm)", ylabel="Number of particles", p2[1])
    plot!(xlabel="Depth (mm)", ylabel="Number of particles", p2[2])
end

title!(p2[1], "N_geom")
title!(p2[2], "N_accept")

display(p2)

# ============================================================
# PLOT 3: FRAZIONE VOLUMICA MONTE CARLO
# ============================================================
pφ = plot(size=(700, 500),
    margin=5Plots.mm,
    left_margin=10Plots.mm,
    bottom_margin=5Plots.mm,
    frame=:box)
colors = [
    :navy,            # blu scuro (riferimento)
    :royalblue,      # blu
    :steelblue,      # blu-grigio
    :lightgray,      # diamond grey
    :deepskyblue,    # azzurro intenso
    :darkorange      # arancio contrasto
]
yticks_positions = 0:5:50
for (i, n) in enumerate(n_arr)
    plot!(pφ, x_mm, all_phi[i] .* 100,
        color=colors[i],
        xticks=xticks_positions,
        yticks=yticks_positions,
        lw=2.5,
        frame=:box,
        gridalpha=0.05,
        label="r=$(n * r_Li * 10000) μm",
        xlim=(0, sat * 10))
end


vline!(pφ, [x_split], ls=:dash, color=:black, label="RDR/2")

xlabel!(pφ, "Depth (mm)")
ylabel!(pφ, "Volume fraction (%)")
display(pφ)
