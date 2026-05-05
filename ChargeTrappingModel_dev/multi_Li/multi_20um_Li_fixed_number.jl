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
# GEOMETRIA
# ============================================================
const Lx = 0.2
const Ly = 0.2
const Lz = 0.2
const dx = 0.001
const FCCD_cm = 0.1
const α = 1.6e-11
# ============================================================
# DIFFUSIONE
# ============================================================
const D = 28.9
const Δt = 1.0
const t_max = 10000
const Nt = Int(t_max / Δt)

const σ = sqrt(6 * D * Δt) * 1e-4

# ============================================================
# LITIO
# ============================================================
const r_Li = 0.002
const cell_size = 2 * 0.012

# ============================================================
# REGIONI
# ============================================================
#const r_regions = (0.012, 0.010, 0.008, 0.006, 0.004, 0.002)
const r_regions = (0.012, 0.011, 0.010, 0.009, 0.008, 0.006, 0.004, 0.002)
const n_regions = length(r_regions)
const region_bounds = collect(range(0, 1, length=n_regions + 1))

# ============================================================
# ANNEALING
# ============================================================
const t_ann = 18 * 60
const T_ann = 623

const R = 1.98
const H = 11800
const D0 = 2.5e-3

const Ns = 10^(21.27 - 2610 / T_ann)
const D_Li = D0 * exp(-H / (R * T_ann))

const Nd_saturation = 1e14
const packing = 0.64




# =====
# FIXED Li
# =====
N_input = [
    3, 3, 3, 3, 3, 3, 3, 3,
    4, 4, 4, 4, 4, 4, 4, 4,
    4, 3, 6, 4, 6, 4, 4, 6,
    7, 8, 5, 8, 7, 7, 4, 3,
    5, 2, 6, 4, 5, 6, 1, 2,
    0, 0, 1, 1, 1, 0, 0, 1,
    1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,]




# ============================================================
# SATURATION DEPTH (precomputed once)
# ============================================================
depth_list = collect(0:dx:0.11)
Nd_vals = Ns .* erfc.(depth_list ./ (2 * sqrt(D_Li * t_ann)))

idx = findfirst(i -> Nd_vals[i] - Nd_saturation <= 0, eachindex(depth_list))
const saturation_depth = depth_list[idx]

println("Saturation depth = ", saturation_depth * 10, " mm")

# ============================================================
# FAST REGION FUNCTION
# ============================================================
@inline function get_region(x)
    ξ = x / saturation_depth
    @inbounds for i in 1:n_regions
        if region_bounds[i] <= ξ <= region_bounds[i+1]
            return i
        end
    end
    return 0
end

# ============================================================
# OVERLAP CHECK (optimized)
# ============================================================
@inline function is_overlapping(x, y, z, r, cells, nx, ny, nz)

    ix = clamp(Int(fld(x, cell_size)) + 1, 1, nx)
    iy = clamp(Int(fld(y, cell_size)) + 1, 1, ny)
    iz = clamp(Int(fld(z, cell_size)) + 1, 1, nz)

    @inbounds for dx in -1:1, dy in -1:1, dz in -1:1
        jx = ix + dx
        jy = iy + dy
        jz = iz + dz

        if 1 ≤ jx ≤ nx && 1 ≤ jy ≤ ny && 1 ≤ jz ≤ nz
            cell = cells[jx, jy, jz]
            @inbounds for p in cell
                dxp = x - p.x
                dyp = y - p.y
                dzp = z - p.z
                if dxp * dxp + dyp * dyp + dzp * dzp < (r + p.r)^2
                    return true
                end
            end
        end
    end
    return false
end

# ============================================================
function truncated_poisson(λ, Nmax)
    if Nmax <= 0
        return 0
    end

    for _ in 1:50
        n = rand(Poisson(λ))
        if n <= Nmax
            return n
        end
    end

    return Nmax
end

# ============================================================
# GENERAZIONE
# ============================================================
function generate_Li()

    nx = floor(Int, Lx / cell_size)
    ny = floor(Int, Ly / cell_size)
    nz = floor(Int, Lz / cell_size)

    cells = [Vector{LiParticle}() for _ in 1:nx, _ in 1:ny, _ in 1:nz]

    x_slices = collect(0:dx:Lx)

    Ni_geom = zeros(Int, length(x_slices))
    Ni_accept = zeros(Int, length(x_slices))



    nx_c = nx
    ny_c = ny
    nz_c = nz

    total = 0

    @inbounds for (i, xi) in enumerate(x_slices)

        if xi >= saturation_depth
            continue
        end

        region = get_region(xi)
        if region == 0
            continue
        end

        r_big = r_regions[region]

        println("Slice $i at x = $(round(xi*10, digits=2)) mm, region = $region, r_big = $(round(r_big*10, digits=3)) μm")
        N_big = N_input[min(i, length(N_input))]

        scale = (r_big / r_Li)^3
        N_i = Int(round(N_big * scale))

        V_particle = (4 / 3) * π * r_Li^3
        V_slice = dx * Ly * Lz
        N_max = Int(floor(packing * V_slice / V_particle))

        N_i = min(N_i, N_max)
        Ni_geom[i] = N_i

        accepted = 0
        tries_max = 25 * N_i

        @inbounds for _ in 1:tries_max
            if accepted >= N_i
                break
            end

            x = xi + rand() * dx
            y = rand() * Ly
            z = rand() * Lz

            if !is_overlapping(x, y, z, r_Li, cells, nx_c, ny_c, nz_c)

                ix = clamp(Int(fld(x, cell_size)) + 1, 1, nx_c)
                iy = clamp(Int(fld(y, cell_size)) + 1, 1, ny_c)
                iz = clamp(Int(fld(z, cell_size)) + 1, 1, nz_c)

                push!(cells[ix, iy, iz], LiParticle(x, y, z, r_Li))

                accepted += 1
                total += 1
                Ni_accept[i] += 1
            end
        end
    end

    println("Total Li particles = $total")

    return cells, x_slices, Ni_accept, Ni_geom
end



# ============================================================
# RUN GENERATION
# ============================================================
cells, x_slices, Ni_accept, Ni_geom = generate_Li()


# ============================================================
# PLOT
# ============================================================
#=
x_mm = x_slices .* 10

p_compare = plot(
    x_mm,
    Ni_geom,
    lw=2,
    color=:black,
    label="Geometric + Poisson limit",
    xlabel="Depth (mm)",
    ylabel="Number of Li per slice",
    frame=:box
)

plot!(
    x_mm,
    Ni_accept,
    lw=2,
    ls=:dash,
    color=:orange,
    label="Accepted Li (simulation)"
)

vline!(
    p_compare,
    [saturation_depth * 10],
    ls=:dash,
    color=:blue,
    label="Saturation depth"
)

display(p_compare)
savefig(p_compare, "plot/Li_geom_vs_accept.png")=#

# ============================================================
# TRANSPORT (OPTIMIZED BUT SAME LOGIC)
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

    σ_local = σ

    @inbounds for n in 1:N

        x = x_charges[n]
        y = rand() * Ly
        z = rand() * Lz

        is_trapped = false

        for _ in 1:Nt

            x += σ_local * randn()
            y += σ_local * randn()
            z += σ_local * randn()

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

                ix = clamp(Int(fld(x, cell_size)) + 1, 1, nx)
                iy = clamp(Int(fld(y, cell_size)) + 1, 1, ny)
                iz = clamp(Int(fld(z, cell_size)) + 1, 1, nz)

                @inbounds for dxi in -1:1, dyi in -1:1, dzi in -1:1
                    jx, jy, jz = ix + dxi, iy + dyi, iz + dzi

                    if 1 ≤ jx ≤ nx && 1 ≤ jy ≤ ny && 1 ≤ jz ≤ nz
                        cell = cells[jx, jy, jz]
                        @inbounds for p in cell
                            dxp = x - p.x
                            dyp = y - p.y
                            dzp = z - p.z
                            if dxp * dxp + dyp * dyp + dzp * dzp < r^2
                                is_trapped = true
                                break
                            end
                        end
                        if is_trapped
                            break
                        end
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
# CCE SIMULATION (UNCHANGED LOGIC)
# ============================================================
x_pos = collect(0:dx:0.11)

N_charges = 250
N_repeat = 25

CCE_matrix = zeros(length(x_pos), N_repeat)

pbar = Progress(length(x_pos) * N_repeat)

@inbounds for (i, x0) in enumerate(x_pos)
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
# STATISTICS
# ============================================================
CCE_mean = vec(mean(CCE_matrix, dims=2))
CCE_std = vec(std(CCE_matrix, dims=2))
xticks_positions = 0:0.1:1.5
yticks_positions = 0:0.05:1
p = plot(
    x_pos .* 10, CCE_mean,
    ribbon=CCE_std,
    lw=2,
    color=:deepskyblue,
    xlabel="depth (mm)",
    ylabel="CCE",
    xticks=(xticks_positions, string.(xticks_positions)),
    yticks=(yticks_positions, string.(yticks_positions)),
    label="CCE",
    frame=:box,
    size=(400, 600)
)

display(p)

# ============================================================
# SAVE JSON
# ============================================================
#filename = "JSON-packing/CCE_8slices_fixed_number.json"
filename = "JSON-packing/CCE_8slices_fixed_number_agglomerate.json"

results = Dict(
    "x_pos" => collect(x_pos),
    "CCE_mean" => CCE_mean,
    "CCE_std" => CCE_std
)

if isfile(filename)
    println("File exists. Overwrite? (y/n)")
    if lowercase(readline()) != "y"
        println("Cancelled.")
        return
    end
end

open(filename, "w") do f
    JSON.print(f, results, 4)
end

println("Saved: $filename")