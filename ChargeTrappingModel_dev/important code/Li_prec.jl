using Random
using Plots
using SpecialFunctions
using Distributions
using JSON
using ProgressMeter

gr()

# ============================================================
# STRUCT
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
# LITIO
# ============================================================
r_Li = 0.002          # 20 μm
cell_size = r_Li      # come richiesto

# ============================================================
# ANNEALING
# ============================================================
t_ann = 18 * 60
T_ann = 623

R = 1.98
H = 11800
D0 = 2.5e-3

Ns = 10^(21.27 - 2610 / T_ann)
D_Li = D0 * exp(-H / (R * T_ann))

Nd_saturation = 1e14
α = 4e-10

# ============================================================
# SATURATION DEPTH
# ============================================================
depth_list = 0:dx:0.11
Nd_vals = Ns .* erfc.(depth_list ./ (2 * sqrt(D_Li * t_ann)))

idx = findfirst(i -> Nd_vals[i] - Nd_saturation <= 0, eachindex(depth_list))
saturation_depth = depth_list[idx]

println("Saturation depth = ", saturation_depth * 10, " mm")

# ============================================================
# POISSON TRONCATA
# ============================================================
function truncated_poisson(λ, Nmax)
    if Nmax <= 0
        return 0
    end

    for _ in 1:100
        n = rand(Poisson(λ))
        if n <= Nmax
            return n
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

    nx = floor(Int, Lx / cell_size)
    ny = floor(Int, Ly / cell_size)
    nz = floor(Int, Lz / cell_size)

    cells = [Vector{LiParticle}() for _ in 1:nx, _ in 1:ny, _ in 1:nz]

    x_slices = collect(0:dx:Lx)

    total = 0

    for xi in x_slices

        if xi >= saturation_depth
            continue
        end

        Nd_val = Ns * erfc(xi / (2 * sqrt(D_Li * t_ann)))
        density = α * max(Nd_val - Nd_saturation, 0.0)

        λ = density * dx * Ly * Lz

        # limite geometrico
        V_particle = (4 / 3) * π * r_Li^3
        V_slice = dx * Ly * Lz
        N_max = Int(floor(V_slice / V_particle))

        N_i = truncated_poisson(λ, N_max)

        trials = N_i * 100
        accepted = 0

        for _ in 1:trials

            if accepted >= N_i
                break
            end

            x = xi + rand() * dx
            y = rand() * Ly
            z = rand() * Lz

            if !is_overlapping(x, y, z, r_Li, cells, nx, ny, nz, cell_size)

                ix = clamp(Int(floor(x / cell_size)) + 1, 1, nx)
                iy = clamp(Int(floor(y / cell_size)) + 1, 1, ny)
                iz = clamp(Int(floor(z / cell_size)) + 1, 1, nz)

                push!(cells[ix, iy, iz], LiParticle(x, y, z, r_Li))

                accepted += 1
                total += 1
            end
        end
    end

    println("Total Li particles = $total")
    return cells
end

# ============================================================
# RUN GENERAZIONE
# ============================================================
cells = generate_Li()

# ============================================================
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
N_charges = 250
N_repeat = 25

CCE_matrix = zeros(length(x_pos), N_repeat)

pbar = Progress(length(x_pos) * N_repeat)

for (i, x0) in enumerate(x_pos)
    for j in 1:N_repeat

        collected, _ =
            multiple_charges_trapping_3D(
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
# STATISTICHE + PLOT
# ============================================================
CCE_mean = vec(mean(CCE_matrix, dims=2))
CCE_std = vec(std(CCE_matrix, dims=2))

p = plot(
    x_pos .* 10, CCE_mean,
    ribbon=CCE_std,
    lw=2,
    xlabel="depth (mm)",
    ylabel="CCE",
    label="CCE",
    frame=:box,
    size=(450, 600)
)

display(p)


filename = "JSON/CCE_20um.json"
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

