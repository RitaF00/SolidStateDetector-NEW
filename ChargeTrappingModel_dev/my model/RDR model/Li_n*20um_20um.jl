using Plots
using SpecialFunctions
using Distributions
using JSON
using ProgressMeter

gr()

# ============================================================
# NOTE FISICHE
# ============================================================
"""
Modello di crescita litio:
- raggio dipende dalla profondità
- riempimento con vincolo geometrico + Poisson
Si inserisceil numero massimo di litio che può essere geenrato tramite
un packing factor stimato essere 0.64 in 3D
[file:///C:/Users/rirri/Downloads/PhysRevA.27.1053.pdf]

In questo codice vi sono anche i plot:
ù- numero massimo di litio permesso geometricamente e quello effettivamente richiesto
- distribuzione 3D dei precipitati
- distribuzione del raggio in funzione della profondità
- plot della frazione del voluumem occupato considerando che alcuni precipitati sono più grandi di altri (quindi non è detto che 1000 Li da 20 μm siano equivalenti a 1 Li da 120 μm)
"""

# ============================================================
# STRUTTURE
# ============================================================
struct LiParticle
    x::Float64
    y::Float64
    z::Float64
    r::Float64
end

# ============================================================
# PARAMETRI GEOMETRICI (cm)
# ============================================================
Lx = 0.2
Ly = 0.2
Lz = 0.2

dx = 0.001
FCCD_cm = 0.1

# ============================================================
# DIFFUSIONE
# ============================================================
D = 28.9
Δt = 1.0
t_max = 10000
Nt = Int(t_max / Δt)

σ = sqrt(6 * D * Δt) * 1e-4

# ============================================================
# TRAPPING
# ============================================================
r_Li = 0.002   # 20 μm in cm

# ============================================================
# PARAMETRO GLOBALE (IMPORTANTE)
# ============================================================
n = 6 # <-- fattore di ingrandimento raggio regione RDR
cell_size = 2 * n * r_Li  #rendiamo cell_size globale
# ============================================================
# ANNEALING / MATERIAL PARAMETERS
# ============================================================
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
depth_list = 0:dx:0.11
Nd_vals = Ns .* erfc.(depth_list ./ (2 * sqrt(D_Li * t_ann)))

idx = findfirst(i -> Nd_vals[i] - Nd_saturation <= 0, eachindex(depth_list))
saturation_depth = depth_list[idx]

println("Saturation depth = ", saturation_depth * 10, " mm")

# ============================================================
# PROFILO RAGGIO
# ============================================================
function r_Li_profile(x)
    return x <= saturation_depth / 2 ? n * r_Li : r_Li
end

# ============================================================
# POISSON TRONCATA
# ============================================================
function truncated_poisson(λ, Nmax)
    if Nmax <= 0
        return 0
    end

    for _ in 1:50
        N = rand(Poisson(λ))
        if N <= Nmax
            return N
        end
    end

    return Nmax
end

# ============================================================
# OVERLAP CHECK
# ============================================================
function is_overlapping(x, y, z, r, cells, nx, ny, nz, cell_size)

    ix = clamp(Int(floor(x / cell_size)) + 1, 1, nx)
    iy = clamp(Int(floor(y / cell_size)) + 1, 1, ny)
    iz = clamp(Int(floor(z / cell_size)) + 1, 1, nz)

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
function generate_Li()

    println(cell_size)

    nx = floor(Int, Lx / cell_size)
    ny = floor(Int, Ly / cell_size)
    nz = floor(Int, Lz / cell_size)


    cells = [Vector{LiParticle}() for _ in 1:nx, _ in 1:ny, _ in 1:nz]

    x_slices = collect(0:dx:Lx)

    Ni_geom = zeros(Int, length(x_slices))
    Ni_accept = zeros(Int, length(x_slices))

    packing = 0.64
    total = 0

    for (i, xi) in enumerate(x_slices)

        if xi >= saturation_depth
            continue
        end

        Nd_val = Ns * erfc(xi / (2 * sqrt(D_Li * t_ann)))
        density = α * max(Nd_val - Nd_saturation, 0.0)

        λ = density * dx * Ly * Lz

        r = r_Li_profile(xi)

        print(r)

        V_particle = (4 / 3) * π * r^3
        V_slice = dx * Ly * Lz
        N_max = Int(floor(packing * V_slice / V_particle))

        N_target = truncated_poisson(λ, N_max)

        Ni_geom[i] = N_target

        accepted = 0
        trials = max(20 * N_target, 50)

        for _ in 1:trials

            if accepted >= N_target
                break
            end

            x = xi + rand() * dx
            y = rand() * Ly
            z = rand() * Lz

            r = r_Li_profile(x)

            if !is_overlapping(x, y, z, r, cells, nx, ny, nz, cell_size)

                ix = clamp(Int(floor(x / cell_size)) + 1, 1, nx)
                iy = clamp(Int(floor(y / cell_size)) + 1, 1, ny)
                iz = clamp(Int(floor(z / cell_size)) + 1, 1, nz)

                push!(cells[ix, iy, iz], LiParticle(x, y, z, r))

                accepted += 1
                total += 1
            end
        end

        Ni_accept[i] = accepted
    end

    println("\nTotal accepted particles = $total")

    return cells, Ni_geom, Ni_accept, x_slices
end


# ============================================================
# RUN
# ============================================================
cells, Ni_geom, Ni_accept, x_slices = generate_Li()

# ============================================================
# STATISTICHE
# ============================================================
println("\n=== STATISTICHE ===")
println("Media richiesti = ", mean(Ni_geom))
println("Media accettati = ", mean(Ni_accept))
println("Efficienza media = ", mean(Ni_accept ./ max.(Ni_geom, 1)))

x_mm = x_slices .* 10
x_split = saturation_depth * 10 / 2
yticks_positions = 0:floor(Int, maximum(Ni_geom) / 5):maximum(Ni_geom)
xticks_positions = 0:0.01:0.7


p = plot()

# richiesti
plot!(p,
    x_mm, Ni_geom,
    xlim=[0, 0.7],
    lw=3,
    yticks=yticks_positions,
    xticks=xticks_positions,
    color=:black,
    gridalpha=0.2,
    label="Requested (Poisson + geom)"
)

# accettati
plot!(p,
    x_mm, Ni_accept,
    lw=2,
    ls=:dash,
    color=:orange,
    label="Accepted (no overlap)"
)

xlabel!("Depth (mm)")
ylabel!("Li precipitates per slice")

vline!(p, [x_split],
    ls=:dash,
    color=:blue,
    label="RDR boundary"
)

plot!(p,
    frame=:box,
    size=(800, 500),
    legend=:topleft,
    legendfontsize=10
)

display(p)


function volume_fraction_MC(cells, x_slices, dx, Lx, Ly, Lz; Nsamp=2000)

    nx, ny, nz = size(cells)

    φ = zeros(length(x_slices))

    V_slice = dx * Ly * Lz

    for (i, xi) in enumerate(x_slices)

        inside = 0

        for _ in 1:Nsamp

            # punto random nella slice
            x = xi + rand() * dx
            y = rand() * Ly
            z = rand() * Lz

            # trova cella
            ix = clamp(Int(floor(x / cell_size)) + 1, 1, nx)
            iy = clamp(Int(floor(y / cell_size)) + 1, 1, ny)
            iz = clamp(Int(floor(z / cell_size)) + 1, 1, nz)

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

                if found
                    break
                end
            end
        end

        φ[i] = inside / Nsamp
    end

    return φ
end
φ = volume_fraction_MC(cells, x_slices, dx, Lx, Ly, Lz, Nsamp=3000)
x_mm = x_slices .* 10
x_split = saturation_depth * 10 / 2

pφ = plot(
    x_mm, φ .* 100,
    xlim=[0, 0.7],
    lw=2,
    color=:purple,
    xlabel="Depth (mm)",
    ylabel="Volume fraction",
    label="φ %",
    frame=:box,
    size=(800, 500)
)

vline!(pφ, [x_split],
    ls=:dash,
    color=:black,
    label="RDR boundary"
)

display(pφ)


# TRASPORTO CARICHE
# ============================================================
function multiple_charges_trapping_3D(
    x_charges, N, cells,
    x_saturation, dx, Lx, Ly, Lz,
    Nt, σ, FCCD_cm
)

    if isa(x_charges, Number)
        x_charges = fill(x_charges, N)
    end

    nx, ny, nz = size(cells)

    collected = 0
    trapped = 0

    for n in 1:N

        x = x_charges[n]
        y = rand() * Ly
        z = rand() * Lz

        is_trapped = false

        #r = r_Li_profile(x)
        #println(r)

        for _ in 1:Nt

            x += σ * randn()
            y += σ * randn()
            z += σ * randn()

            x = clamp(x, 0, Lx)
            y = clamp(y, 0, Ly)
            z = clamp(z, 0, Lz)

            if x >= FCCD_cm
                collected += 1
                break
            end

            if x <= dx
                is_trapped = true
                break
            end

            if x <= x_saturation

                ix = clamp(Int(floor(x / cell_size)) + 1, 1, nx)
                iy = clamp(Int(floor(y / cell_size)) + 1, 1, ny)
                iz = clamp(Int(floor(z / cell_size)) + 1, 1, nz)

                for dxi in -1:1, dyi in -1:1, dzi in -1:1

                    jx, jy, jz = ix + dxi, iy + dyi, iz + dzi

                    if 1 ≤ jx ≤ nx && 1 ≤ jy ≤ ny && 1 ≤ jz ≤ nz

                        for p in cells[jx, jy, jz]
                            if (x - p.x)^2 + (y - p.y)^2 + (z - p.z)^2 < p.r^2
                                is_trapped = true
                                break
                            end
                        end
                    end

                    if is_trapped
                        break
                    end
                end
            end

            if is_trapped
                trapped += 1
                break
            end
        end
    end

    return collected, trapped
end

# ============================================================
# SIMULAZIONE CCE
# ============================================================
x_pos = 0:dx:0.11
N_charges = 150
N_repeat = 15




CCE_matrix = zeros(length(x_pos), N_repeat)

pbar = Progress(length(x_pos) * N_repeat)

for (i, x0) in enumerate(x_pos)
    for j in 1:N_repeat

        collected, trapped = multiple_charges_trapping_3D(
            x0, N_charges, cells,
            saturation_depth, dx,
            Lx, Ly, Lz,
            Nt, σ,
            FCCD_cm
        )

        CCE_matrix[i, j] = collected / N_charges
        next!(pbar)
    end
end

# ============================================================
# STATISTICS
# ============================================================
CCE_mean = vec(mean(CCE_matrix, dims=2))
CCE_std = vec(std(CCE_matrix, dims=2))

# ============================================================
# PLOT (PAPER STYLE)
# ============================================================
xticks_positions = 0:0.1:1.1

p = scatter(
    x_pos .* 10,
    CCE_mean,
    yerror=CCE_std,
    markersize=3,
    label="CCE",
    xlabel="Depth (mm)",
    ylabel="CCE",
    frame=:box,
    xticks=(xticks_positions, string.(round.(xticks_positions, digits=2))),
    ylim=(0, 1),
    size=(500, 400),
    dpi=300
)

display(p)

# ============================================================
# JSON SAVE
# ============================================================
filename = "JSON-packing/CCE_$(n*r_Li*1e4)um_20um.json"
#filename = "JSON/CCE_prova.json"
results = Dict(
    "x_pos" => collect(x_pos),
    "CCE_mean" => CCE_mean,
    "CCE_std" => CCE_std,
    "n_factor" => n,
    "r_Li_cm" => r_Li,
    "saturation_depth" => saturation_depth
)

mkpath(dirname(filename))

if isfile(filename)
    println("File esiste. Sovrascrivere? (y/n)")
    ans = readline()
    if ans != "y" && ans != "Y"
        println("Abort.")
        return
    end
end

open(filename, "w") do f
    JSON.print(f, results, 4)
end

println("Salvato: $filename")

#=


# ============================================================
# RADIUS vs DEPTH PLOT
# ============================================================
x = collect(0:dx:Lx)

r_vals = [r_Li_profile(xi) for xi in x]

xticks_positions = 0:0.1:1.1

x = collect(0:dx:saturation_depth)

r_vals = [r_Li_profile(xi) for xi in x]
cell_vals = [2 * n * r_Li for xi in x]

xticks_positions = 0:0.1:1.1

p_r = plot(
    x .* 10,
    r_vals .* 1e4,   # cm → μm
    lw=2,
    color=:black,
    xlabel="Depth (mm)",
    ylabel="Size (μm)",
    label="r_Li(x)",
    frame=:box,
    size=(600, 400),
    xticks=(xticks_positions, string.(round.(xticks_positions, digits=2))),
    grid=false
)
#=
plot!(
    p_r,
    x .* 10,
    cell_vals .* 1e4,
    lw=2,
    color=:blue,
    label="cell_size(x)"
)=#

# transizione RDR
vline!(p_r, [saturation_depth * 10 / 2],
    ls=:dash,
    color=:red,
    label="RDR boundary"
)








# =========================
# EXTRACT POINTS
# =========================
xs, ys, zs, sizes, colors = Float64[], Float64[], Float64[], Float64[], Any[]

for cell in cells
    for p in cell
        push!(xs, p.x)
        push!(ys, p.y)
        push!(zs, p.z)

        if p.x <= saturation_depth / 2
            push!(sizes, 6)
            push!(colors, :blue)
        else
            push!(sizes, 2)
            push!(colors, :red)
        end
    end
end

# =========================
# 3D PLOT
# =========================
p3d = scatter3d(xs, ys, zs,
    markersize=sizes,
    markercolor=colors,
    alpha=0.7,
    xlabel="x",
    ylabel="y",
    zlabel="z",
    legend=false,
    size=(800, 800),
    dpi=300
)
xv = [0, Lx]
yv = [0, Ly]
zv = [0, Lz]

# Aggiungo le 12 linee degli spigoli
for x in xv, y in yv
    plot!(p3d, [x, x], [y, y], [0, Lz], linecolor=:black)  # spigoli verticali
end
for x in xv, z in zv
    plot!(p3d, [x, x], [0, Ly], [z, z], linecolor=:black)  # spigoli front/back
end
for y in yv, z in zv
    plot!(p3d, [0, Lx], [y, y], [z, z], linecolor=:black)  # spigoli sinistra/destra
end
display(p3d)
savefig(p3d, "li_120um_20um/ 3D_dist.png")

# =========================
# PROJECTIONS
# =========================
p_xy = scatter(xs, ys,
    markersize=sizes,
    markercolor=colors,
    xlabel="x",
    ylabel="y",
    legend=false,
    title="XY",
    frame=:box
)


p_xz = scatter(xs, zs,
    markersize=sizes,
    markercolor=colors,
    xlabel="x",
    ylabel="z",
    legend=false,
    title="XZ",
    frame=:box
)

p_yz = scatter(ys, zs,
    markersize=sizes,
    markercolor=colors,
    xlabel="y",
    ylabel="z",
    legend=false,
    title="YZ",
    frame=:box
)
p = plot(p_xy, p_xz, p_yz, layout=(1, 3), dpi=300, size=(2000, 600))
savefig(p, "li_120um_20um/ 3D_dist_xy_xz_zy.png")


# =========================
# SLICE PLOT (IMPORTANT FIX)
# =========================

x_mm = x_slices .* 10
x_split = saturation_depth * 10 / 2

idx1 = x_mm .<= x_split
idx2 = x_mm .>= x_split

p = plot()

# prima parte (blu)
plot!(p,
    x_mm[idx1], Ni_geom[idx1],
    lw=2,
    color=:deepskyblue,
    label="120 μm",
    xlim=[0, 0.7]
)

# seconda parte (arancione)
plot!(p,
    x_mm[idx2], Ni_geom[idx2],
    lw=2,
    color=:orange,
    label="20 μm"
)

xlabel!("x (mm)")
ylabel!("Li per slice")

plot!(p,
    legend=:topright,
    frame=:box,
    size=(700, 500)
)

vline!(p, [x_split],
    ls=:dash,
    color=:grey,
    alpha=0.3,
    label="RDR/2"
)
savefig(p, "li_120um_20um/Li_sites_each_slice.png")=#


