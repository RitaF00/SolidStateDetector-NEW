using Random
using Plots
using SpecialFunctions
using Distributions
using JSON
using ProgressMeter

gr()

# =====================
# FUNZIONE r_Li(x)
# =====================
function r_Li_depth(x, x_saturation)
    if x <= x_saturation / 2
        return 0.02   # 200 μm in cm
    elseif x <= x_saturation
        return 0.002  # 20 μm in cm
    else
        return 0.0
    end
end


# =====================
# Li generation
# =====================
function generate_Li_cells(Lx, Ly, Lz, dx, α, Ns, D_Li, t_ann, x_saturation)

    cell_size = 0.002  # 20 μm

    nx = Int(Lx / cell_size)
    ny = Int(Ly / cell_size)
    nz = Int(Lz / cell_size)

    cells = [Vector{NTuple{3,Float64}}() for _ in 1:nx, _ in 1:ny, _ in 1:nz]

    Nd_saturation = 1e14
    x_slices = 0:dx:Lx
    total_Li = 0

    for xi in x_slices

        if xi >= x_saturation
            continue
        end

        Nd_val = Ns * erfc(xi / (2 * sqrt(D_Li * t_ann)))
        effective_Nd = max(Nd_val - Nd_saturation, 0.0)

        dV = dx * Ly * Lz
        λ = α * effective_Nd * dV

        N_i = rand(Poisson(λ))
        total_Li += N_i

        for _ in 1:N_i
            x_pos = xi + rand() * dx
            y_pos = rand() * Ly
            z_pos = rand() * Lz

            ix = clamp(Int(floor(x_pos / cell_size)) + 1, 1, nx)
            iy = clamp(Int(floor(y_pos / cell_size)) + 1, 1, ny)
            iz = clamp(Int(floor(z_pos / cell_size)) + 1, 1, nz)

            push!(cells[ix, iy, iz], (x_pos, y_pos, z_pos))
        end
    end

    println("Total precipitates = ", total_Li)
    return cells, nx, ny, nz, cell_size
end


# =====================
# Charge transport + trapping
# =====================
function multiple_charges_trapping_3D(
    x_charges, N, cells, cell_size,
    x_saturation, dx, Lx, Ly, Lz,
    Nt, σ, FCCD_cm
)

    if isa(x_charges, Number)
        x_charges = fill(x_charges, N)
    end

    nx, ny, nz = size(cells)

    collected_count = 0
    trapped_count = 0

    # 🔥 search radius per precipitati grandi (200 μm)
    max_r = 0.02
    search_radius = ceil(Int, max_r / cell_size)

    for n in 1:N

        x = x_charges[n]
        y = rand() * Ly
        z = rand() * Lz

        trapped = false

        for _ in 1:Nt

            # diffusion
            x += σ * randn()
            y += σ * randn()
            z += σ * randn()

            # boundaries
            x = clamp(x, 0, Lx)
            y = clamp(y, 0, Ly)
            z = clamp(z, 0, Lz)

            # collection
            if x >= FCCD_cm
                collected_count += 1
                break
            end

            # loss at entrance
            if x <= dx
                trapped = true
                break
            end

            # trapping region
            if x <= x_saturation

                ix = clamp(Int(floor(x / cell_size)) + 1, 1, nx)
                iy = clamp(Int(floor(y / cell_size)) + 1, 1, ny)
                iz = clamp(Int(floor(z / cell_size)) + 1, 1, nz)

                for dxi in -search_radius:search_radius,
                    dyi in -search_radius:search_radius,
                    dzi in -search_radius:search_radius

                    jx, jy, jz = ix + dxi, iy + dyi, iz + dzi

                    if 1 ≤ jx ≤ nx && 1 ≤ jy ≤ ny && 1 ≤ jz ≤ nz

                        for (px, py, pz) in cells[jx, jy, jz]

                            r_local = r_Li_depth(px, x_saturation)

                            if (x - px)^2 + (y - py)^2 + (z - pz)^2 < r_local^2
                                trapped = true
                                break
                            end
                        end
                    end

                    if trapped
                        break
                    end
                end
            end

            if trapped
                trapped_count += 1
                break
            end
        end
    end

    return collected_count, trapped_count
end


# =====================
# PARAMETERS
# =====================
Lx = 0.2
Ly = 0.2
Lz = 0.2

dx = 0.001
FCCD_cm = 0.1

# diffusion
D = 28.9
Δt = 1.0
t_max = 10000
Nt = Int(t_max / Δt)
σ = sqrt(6 * D * Δt) * 1e-4

# annealing
t_ann = 18 * 60
T_ann = 623

R = 1.98
H = 11800
D0 = 2.5e-3

Ns = 10^(21.27 - 2610 / T_ann)
D_Li = D0 * exp(-H / (R * T_ann))

Nd_saturation = 1e14

# =====================
# Saturation depth
# =====================
depth_list = 0:dx:0.11
Nd(x) = Ns .* erfc.(x ./ (2 * sqrt(D_Li * t_ann)))
Nd_vals = Nd(depth_list)

diff_vals = Nd_vals .- Nd_saturation
idx = findfirst(i -> diff_vals[i] <= 0, eachindex(diff_vals))
saturation_depth = depth_list[idx]

println("=============================================================================")
println("Saturation depth = ", saturation_depth * 10, " mm")
println("=============================================================================")


# =====================
# GENERATE DEFECTS
# =====================
α = 1.6e-11

cells, nx, ny, nz, cell_size = generate_Li_cells(
    Lx, Ly, Lz, dx, α, Ns, D_Li, t_ann, saturation_depth
)

println("=============================================================================")


# =====================
# CCE SIMULATION
# =====================
x_pos = 0:dx:0.11
N_charges = 250
N_repeat = 10

CCE_matrix = zeros(length(x_pos), N_repeat)

pbar = Progress(length(x_pos) * N_repeat)

for (i, x0) in enumerate(x_pos)
    for j in 1:N_repeat

        collected, trapped =
            multiple_charges_trapping_3D(
                x0, N_charges, cells, cell_size,
                saturation_depth, dx,
                Lx, Ly, Lz,
                Nt, σ,
                FCCD_cm
            )

        CCE_matrix[i, j] = collected / N_charges
        next!(pbar)
    end
end


# =====================
# STATISTICS
# =====================
CCE_mean = vec(mean(CCE_matrix, dims=2))
CCE_std = vec(std(CCE_matrix, dims=2))


# =====================
# PLOT
# =====================
p = plot(x_pos .* 10, CCE_mean,
    ribbon=CCE_std,
    lw=2,
    xlabel="depth (mm)",
    ylabel="CCE",
    label="CCE",
    frame=:box,
    size=(450, 600)
)


# =====================
# SAVE JSON
# =====================
filename = "JSON/Nd(x)-Nsat_variable_radius.json"

results = Dict(
    "x_pos" => collect(x_pos),
    "CCE_mean" => CCE_mean,
    "CCE_std" => CCE_std
)

if isfile(filename)
    println("Il file $filename esiste già. Vuoi sovrascriverlo? (y/n)")
    risposta = readline()

    if lowercase(risposta) != "y"
        println("Operazione annullata.")
        return
    end
end

open(filename, "w") do f
    JSON.print(f, results, 4)
end

println("File JSON salvato: $filename")

display(p)
savefig(p, "CCE_variable_radius.png")