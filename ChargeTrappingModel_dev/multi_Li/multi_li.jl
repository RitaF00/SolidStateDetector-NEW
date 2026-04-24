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
In questo caso volgio idvidere l'RDR in più regioni con diversa dimensione di litio
in modo da vedere qual è l'effetto di un gradiente di dimensione del litio (e quindi di diffusività) sulla CCE.
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
n = 6   # <-- fattore di ingrandimento raggio regione RDR
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
#=
function r_Li_profile(x)
    if x <= saturation_depth / 6
        return 0.012
    elseif x <= saturation_depth / 3
        return 0.012
    elseif x <= saturation_depth / 2
        return 0.012
    elseif x <= saturation_depth
        return 0.012
    else
        return 0.002
    end
end
=#

function r_Li_profile(x)
    if x <= saturation_depth / 2
        return r_Li
    else
        return r_Li
    end
end

# ============================================================
# POISSON TRONCATA
# ============================================================
"""
In questo caso, per evitare di generare un numero di particelle troppo grande in alcune slice,
 applichiamo una Poisson troncata che limita il numero massimo di particelle a Nmax.
 dato lambda e Nmax, la funzione prova a generare un numero da Poisson(lambda) fino a 100 volte.
Se il numero generato è minore o uguale a Nmax, lo restituisce. Altrimenti, dopo 100 tentativi, restituisce Nmax.
"""
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
    r_Li = 0.002
    cell_size = 2 * n * r_Li
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

        # è il centro dell apoisooniana che poi alcolo in truncated_poisson
        λ = density * dx * Ly * Lz

        r = r_Li_profile(xi)

        V_particle = (4 / 3) * π * r^3
        V_slice = dx * Ly * Lz
        N_max = Int(floor(V_slice / V_particle))

        N_target = truncated_poisson(λ, N_max)

        Ni_geom[i] = N_target

        accepted = 0
        trials = N_target * 100

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
    Nt, σ, FCCD_cm, r
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
                            if (x - p.x)^2 + (y - p.y)^2 + (z - p.z)^2 < r^2
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
N_charges = 50
N_repeat = 10

CCE_matrix = zeros(length(x_pos), N_repeat)

pbar = Progress(length(x_pos) * N_repeat)

for (i, x0) in enumerate(x_pos)
    for j in 1:N_repeat

        collected, trapped = multiple_charges_trapping_3D(
            x0, N_charges, cells,
            saturation_depth, dx,
            Lx, Ly, Lz,
            Nt, σ,
            FCCD_cm, r_Li
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
filename = "JSON/CCE_4_slice.json"

results = Dict(
    "x_pos" => collect(x_pos),
    "CCE_mean" => CCE_mean,
    "CCE_std" => CCE_std,
    "n_factor" => n,
    "r_Li_cm" => r_Li,
    "saturation_depth" => saturation_depth
)

if isfile(filename)
    println("File esiste. Sovrascrivere? (y/n)")
    ans = readline()
    if lowercase(ans) != "y"
        println("Abort.")
        return
    end
end

open(filename, "w") do f
    JSON.print(f, results, 4)
end

println("Salvato: $filename")




#=
# fattore di scala (da tarare in base al tuo dominio)
k_size = 1500.0

xs = Float64[]
ys = Float64[]
zs = Float64[]
sizes = Float64[]
colors = RGB[]

for cell in cells
    for p in cell
        push!(xs, p.x * 10)
        push!(ys, p.y * 10)
        push!(zs, p.z * 10)

        # dimensione proporzionale al raggio fisico
        push!(sizes, k_size * p.r)

        # colore sempre legato al raggio normalizzato (ok lasciarlo così)
        r_norm = p.r / (n * r_Li)
        push!(colors, RGB(r_norm, 0.2, 1 - r_norm))
    end
end

display(p3d)






p_xy = scatter(xs, ys,
    markersize=sizes,
    markercolor=colors,
    alpha=0.8,
    legend=false,
    xlabel="x (mm)",
    ylabel="y (mm)",
    title="XY plane"
)

p_xz = scatter(xs, zs,
    markersize=sizes,
    markercolor=colors,
    alpha=0.8,
    legend=false,
    xlabel="x (mm)",
    ylabel="z (mm)",
    title="XZ plane"
)

p_zy = scatter(ys, zs,
    markersize=sizes,
    markercolor=colors,
    alpha=0.8,
    legend=false,
    xlabel="y (mm)",
    ylabel="z (mm)",
    title="YZ plane"
)

plot(p_xy, p_xz, p_zy, layout=(1, 3), size=(1200, 400))=#