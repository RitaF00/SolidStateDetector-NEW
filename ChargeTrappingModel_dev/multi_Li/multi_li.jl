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
#n = 1.5  # <-- fattore di ingrandimento raggio regione RDR
cell_size = 0.024
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

    if x <= saturation_depth / 6
        return 0.012

    elseif x <= 2 * saturation_depth / 6
        return 0.010

    elseif x <= 3 * saturation_depth / 6
        return 0.008

    elseif x <= 4 * saturation_depth / 6
        return 0.006

    elseif x <= 5 * saturation_depth / 6
        return 0.004

    elseif x <= saturation_depth
        return 0.002

    else
        return 0.0  # fuori dominio
    end
end



# ============================================================
# POISSON TRONCATA
# ============================================================
function truncated_poisson(λ, Nmax)
    if Nmax <= 0
        return 0
    end

    for _ in 1:100
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
        N_max = Int(floor(V_slice / V_particle))

        N_target = truncated_poisson(λ, N_max)

        Ni_geom[i] = N_target

        accepted = 0
        trials = N_target * 3

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
N_charges = 200
N_repeat = 20
#=

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
filename = "JSON/CCE_6slices.json"
#filename = "JSON/CCE_prova.json"
results = Dict(
    "x_pos" => collect(x_pos),
    "CCE_mean" => CCE_mean,
    "CCE_std" => CCE_std,
    "r_Li_cm" => r_Li,
    "saturation_depth" => saturation_depth
)

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
=#



# ============================================================
# INTERVALLI GENERALIZZATI
# ============================================================
intervals = [
    (0 / 6, 1 / 6, 0.012, :deepskyblue, 16),
    (1 / 6, 2 / 6, 0.010, :slateblue, 10 / 0.65),
    (2 / 6, 3 / 6, 0.008, :indigo, 8 / 0.65),
    (3 / 6, 4 / 6, 0.006, :purple, 6 / 0.65),
    (4 / 6, 5 / 6, 0.004, :orange, 4 / 0.65),
    (5 / 6, 6 / 6, 0.002, :red, 2 / 0.65)
]

# ============================================================
# FUNZIONE PROPRIETÀ
# ============================================================
function get_properties(x, saturation_depth, intervals)
    for (a, b, r, color, size) in intervals
        if a * saturation_depth <= x <= b * saturation_depth
            return r, color, size
        end
    end
    return 0.0, :black, 1.0
end

# ============================================================
# ESTRAZIONE DATI
# ============================================================
xs, ys, zs, sizes, colors = Float64[], Float64[], Float64[], Float64[], Any[]

for cell in cells
    for p in cell
        push!(xs, p.x)
        push!(ys, p.y)
        push!(zs, p.z)

        _, c, s = get_properties(p.x, saturation_depth, intervals)

        push!(sizes, s)
        push!(colors, c)
    end
end

# ============================================================
# 3D PLOT
# ============================================================
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

# box
xv = [0, Lx]
yv = [0, Ly]
zv = [0, Lz]

for x in xv, y in yv
    plot!(p3d, [x, x], [y, y], [0, Lz], linecolor=:black)
end
for x in xv, z in zv
    plot!(p3d, [x, x], [0, Ly], [z, z], linecolor=:black)
end
for y in yv, z in zv
    plot!(p3d, [0, Lx], [y, y], [z, z], linecolor=:black)
end

display(p3d)

# ============================================================
# PROIEZIONI 2D
# ============================================================

# XY
p_xy = scatter(
    xs, ys, frame=:box,
    markersize=sizes,
    markercolor=colors,
    alpha=0.6,
    xlabel="x",
    ylabel="y",
    title="XY projection",
    legend=false,
    size=(600, 600)
)


# XZ
p_xz = scatter(
    xs, zs,
    markersize=sizes, frame=:box,
    markercolor=colors,
    alpha=0.6,
    xlabel="x",
    ylabel="z",
    title="XZ projection",
    legend=false,
    size=(600, 600)
)


# YZ
p_yz = scatter(
    ys, zs,
    markersize=sizes,
    markercolor=colors, frame=:box,
    alpha=0.6,
    xlabel="y",
    ylabel="z",
    title="YZ projection",
    legend=false,
    size=(600, 600)
)
p = plot(p_xy, p_xz, p_yz, layout=(1, 3), dpi=300, size=(2000, 600))
display(p)

# ============================================================
# RAGGIO vs PROFONDITÀ
# ============================================================
x = collect(0:dx:saturation_depth)

r_vals = [get_properties(xi, saturation_depth, intervals)[1] for xi in x]

p_r = plot(
    x .* 10,
    r_vals .* 1e4,
    lw=2,
    color=:black,
    xlabel="Depth (mm)",
    ylabel="Size (μm)",
    label="r_Li(x)",
    frame=:box,
    size=(600, 400),
    dpi=300
)

# linee intervalli
for (a, _, _, _, _) in intervals[2:end]
    vline!(p_r,
        [a * saturation_depth * 10],
        ls=:dash,
        color=:grey,
        label=""
    )
end

display(p_r)

# ============================================================
# Ni PER SLICE
# ============================================================
x_mm = x_slices .* 10

p2 = plot()

for (a, b, _, color, _) in intervals
    idx = (x_mm .>= a * saturation_depth * 10) .&
          (x_mm .<= b * saturation_depth * 10)

    plot!(p2,
        x_mm[idx],
        Ni_geom[idx],
        lw=2,
        color=color, frame=:box,
        label=""
    )
end

xlabel!("x (mm)")
ylabel!("Li per slice")

vline!(p2,
    [saturation_depth * 10 / 2],
    ls=:dash,
    color=:grey,
    alpha=0.3,
    label="RDR/2"
)

display(p2)