using Plots
using SpecialFunctions
using Distributions
using JSON
using ProgressMeter

gr()


"""
In questo codice guardo quanti precipitati di Li posso mettere in una fetta di 10 μm facendoli variare di raggio.
fino a RDR/2 riemoio con raggi molto grnadi mentre dopo con raggi più piccoli. In questo modo vedo qual è il numero massimo di precipitati che posso mettere
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
# PARAMETRI GLOBALI
# ============================================================
const Lx = 0.2
const Ly = 0.2
const Lz = 0.2

const dx = 0.001
const FCCD_cm = 0.1

# Diffusione
const D = 28.9
const Δt = 1.0
const t_max = 10000
const Nt = Int(t_max / Δt)
const σ = sqrt(6 * D * Δt) * 1e-4

# Litio
const r_Li = 0.002   # cm

# Scaling raggio RDR
n_arr = [1, 1.5, 2, 3, 4, 5, 6]

# Annealing
const t_ann = 18 * 60
const T_ann = 623

const R = 1.98
const H = 11800
const D0 = 2.5e-3

const Ns = 10^(21.27 - 2610 / T_ann)
const D_Li = D0 * exp(-H / (R * T_ann))

const Nd_saturation = 1e14
const α = 1.6e-11

# ============================================================
# PROFILI
# ============================================================
function compute_saturation_depth()
    depth_list = 0:dx:0.11
    Nd_vals = Ns .* erfc.(depth_list ./ (2 * sqrt(D_Li * t_ann)))

    idx = findfirst(i -> Nd_vals[i] - Nd_saturation <= 0, eachindex(depth_list))
    return depth_list[idx]
end

function r_Li_profile(x, n, saturation_depth)
    return x <= saturation_depth / 2 ? n * r_Li : r_Li
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
function generate_Li(n, saturation_depth)

    cell_size = 2 * n * r_Li
    println("Cell size = ", cell_size)

    nx = floor(Int, Lx / cell_size)
    ny = floor(Int, Ly / cell_size)
    nz = floor(Int, Lz / cell_size)

    cells = [Vector{LiParticle}() for _ in 1:nx, _ in 1:ny, _ in 1:nz]

    x_slices = collect(0:dx:Lx)

    Ni_geom = zeros(Int, length(x_slices))
    Ni_accept = zeros(Int, length(x_slices))

    total = 0
    total_rejected = 0

    for (i, xi) in enumerate(x_slices)

        if xi >= saturation_depth
            continue
        end

        Nd_val = Ns * erfc(xi / (2 * sqrt(D_Li * t_ann)))
        density = α * max(Nd_val - Nd_saturation, 0.0)

        λ_neat = density * dx * Ly * Lz

        r = r_Li_profile(xi, n, saturation_depth)

        V_particle = (4 / 3) * π * r^3
        V_slice = dx * Ly * Lz
        N_max = Int(floor(V_slice / V_particle))

        λ = rand(Poisson(λ_neat))
        N_target = min(λ, N_max)

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

            r = r_Li_profile(x, n, saturation_depth)

            if !is_overlapping(x, y, z, r, cells, nx, ny, nz, cell_size)

                ix = clamp(Int(floor(x / cell_size)) + 1, 1, nx)
                iy = clamp(Int(floor(y / cell_size)) + 1, 1, ny)
                iz = clamp(Int(floor(z / cell_size)) + 1, 1, nz)

                push!(cells[ix, iy, iz], LiParticle(x, y, z, r))

                accepted += 1
                total += 1
            else
                total_rejected += 1
            end
        end

        Ni_accept[i] = accepted
    end

    return cells, Ni_geom, Ni_accept, x_slices, total
end

# ============================================================
# PLOTTING
# ============================================================
function run_simulation()

    saturation_depth = compute_saturation_depth()
    println("Saturation depth = ", saturation_depth * 10, " mm")

    p1 = plot()
    p_occ = plot()

    for (i, n) in enumerate(n_arr)

        cells, Ni_sampled, Ni_accept, x_slices, total =
            generate_Li(n, saturation_depth)

        x_mm = x_slices .* 10
        x_split = saturation_depth * 10 / 2

        idx1 = x_mm .<= x_split
        idx2 = x_mm .>= x_split

        plot!(p1,
            x_mm[idx1], Ni_sampled[idx1],
            lw=2,
            label="$(n*r_Li*1e4) μm, Ntot = $total",
            xlim=[0, 0.7],
            legendsize=25
        )

        plot!(p1,
            x_mm[idx2], Ni_sampled[idx2],
            lw=2,
            color=:orange,
            label=""
        )

        rmax = r_Li_profile.(x_slices, n, saturation_depth)

        V_particle = (4 / 3) .* π .* rmax .^ 3
        V_slice = dx * Ly * Lz
        #=
        rmax_single = n * r_Li
        V_particle_single = (4 / 3) * π * rmax_single^3

        f_occ = (Ni_sampled .* V_particle_single) ./ V_slice
        =#

        r_local = r_Li_profile.(x_slices, n, saturation_depth)
        V_particle_local = (4 / 3) .* π .* r_local .^ 3

        f_occ = (Ni_sampled .* V_particle_local) ./ V_slice

        plot!(p_occ,
            x_mm,
            f_occ .* 100,
            lw=2,
            label="$(Int(n*r_Li*1e4)) μm",
            xlim=(0, 0.7),)


    end

    vline!(p1, [0.65 / 2], ls=:dash, color=:grey, alpha=0.3, label="RDR/2")

    xlabel!(p1, "x (mm)")
    ylabel!(p1, "Li per slice")

    plot!(p1, legend=:topright, legendfontsize=15, frame=:box, size=(700, 500))
    display(p1)

    vline!(p_occ, [0.65 / 2], ls=:dash, color=:grey, alpha=0.3, label="RDR/2")

    xlabel!(p_occ, "x (mm)")
    ylabel!(p_occ, "Fraction of occupied volume [%]")

    plot!(p_occ, legend=:topright, frame=:box, legendfontsize=15, size=(700, 500))
    display(p_occ)

end

# ============================================================
# MAIN ENTRY
# ============================================================
function main()
    run_simulation()
end

main()