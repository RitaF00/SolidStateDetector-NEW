using Random
using Plots
using SpecialFunctions
using Distributions
using JSON
using ProgressMeter

gr()

# =========================
# STRUCT
# =========================
struct LiParticle
    x::Float64
    y::Float64
    z::Float64
    r::Float64
end

# =========================
# RADIUS PROFILE
# =========================
function r_Li_profile(x)
    r_max = 20e-4   # 120 µm
    r_min = 20e-4    # 20 µm

    x0 = 0.03  # 0.3 mm in cm

    if x <= x0
        return r_max
    else
        return r_min
    end
end

# =========================
# TRUNCATED POISSON
# =========================
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

# =========================
# GENERATION
# =========================
function generate_Li_cells(
    Lx, Ly, Lz,
    dx, α,
    Ns, D_Li, t_ann,
    x_saturation, r
)

    cell_size = 0.005  # cm (50 µm)

    nx = Int(Lx / cell_size)
    ny = Int(Ly / cell_size)
    nz = Int(Lz / cell_size)

    cells = [Vector{LiParticle}() for _ in 1:nx, _ in 1:ny, _ in 1:nz]

    Nd_saturation = 1e14

    x_slices = 0:dx:Lx

    Ni_arr = Float64[]

    total_Li = 0

    for xi in x_slices

        if xi >= x_saturation
            continue
        end

        # =========================
        # ND PROFILE
        # =========================
        Nd_val = Ns * erfc(xi / (2 * sqrt(D_Li * t_ann)))
        effective_Nd = max(Nd_val - Nd_saturation, 0.0)

        dV = dx * Ly * Lz

        λ = α * effective_Nd * dV

        # =========================
        # PARTICLE SIZE
        # =========================
        #r = r_Li_profile(xi)
        V_li = (4 / 3) * π * r^3

        V_slice = dV

        f = 0.1

        N_max = floor(Int, f * V_slice / V_li)
        #println(N_max)

        # Poisson constrained
        N_i = truncated_poisson(λ, N_max)

        push!(Ni_arr, N_i)

        total_Li += N_i

        # =========================
        # PARTICLE GENERATION
        # =========================
        for _ in 1:N_i

            x_pos = xi + rand() * dx
            y_pos = rand() * Ly
            z_pos = rand() * Lz

            p = LiParticle(x_pos, y_pos, z_pos, r)

            ix = clamp(Int(floor(x_pos / cell_size)) + 1, 1, nx)
            iy = clamp(Int(floor(y_pos / cell_size)) + 1, 1, ny)
            iz = clamp(Int(floor(z_pos / cell_size)) + 1, 1, nz)

            push!(cells[ix, iy, iz], p)
        end
    end

    println("Total precipitates = ", total_Li)

    return cells, nx, ny, nz, cell_size, Ni_arr
end

# =========================
# CCE SIMULATION
# =========================
function multiple_charges_trapping_3D(
    x_charges, N,
    cells, cell_size,
    x_saturation, dx,
    Lx, Ly, Lz,
    Nt, σ,
    FCCD_cm
)

    if isa(x_charges, Number)
        x_charges = fill(x_charges, N)
    end

    nx, ny, nz = size(cells)

    collected_count = 0
    trapped_count = 0

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

            # entrance loss
            if x <= dx
                trapped = true
                break
            end

            # trapping
            if x <= x_saturation

                ix = clamp(Int(floor(x / cell_size)) + 1, 1, nx)
                iy = clamp(Int(floor(y / cell_size)) + 1, 1, ny)
                iz = clamp(Int(floor(z / cell_size)) + 1, 1, nz)

                for dxi in -1:1, dyi in -1:1, dzi in -1:1

                    jx, jy, jz = ix + dxi, iy + dyi, iz + dzi

                    if 1 ≤ jx ≤ nx && 1 ≤ jy ≤ ny && 1 ≤ jz ≤ nz

                        for p in cells[jx, jy, jz]

                            if (x - p.x)^2 + (y - p.y)^2 + (z - p.z)^2 < p.r^2
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

# =========================
# PARAMETERS
# =========================
Lx = 0.2
Ly = 0.2
Lz = 0.2

dx = 0.001

FCCD_cm = 0.1

D = 28.9
Δt = 1.0
t_max = 10000
Nt = Int(t_max / Δt)
σ = sqrt(6 * D * Δt) * 1e-4

t_ann = 18 * 60
T_ann = 623

R = 1.98
H = 11800
D0 = 2.5e-3

Ns = 10^(21.27 - 2610 / T_ann)
D_Li = D0 * exp(-H / (R * T_ann))

Nd_saturation = 1e14

depth_list = 0:dx:0.11

Nd(x) = Ns .* erfc.(x ./ (2 * sqrt(D_Li * t_ann)))
Nd_vals = Nd(depth_list)

diff_vals = Nd_vals .- Nd_saturation

idx = findfirst(i -> diff_vals[i] <= 0, eachindex(diff_vals))
saturation_depth = depth_list[idx]

println("Saturation depth = ", saturation_depth * 10, " mm")

# =========================
# GENERATE DEFECTS
# =========================
α = 1.6e-11

r_list = [0.002, 0.006, 0.012]

Ni_all = []
labels = []
for r in r_list
    println("r = $(r*1e4) μm")
    global cells, nx, ny, nz, cell_size, Ni_arr
    cells, nx, ny, nz, cell_size, Ni_arr =
        generate_Li_cells(Lx, Ly, Lz, dx, α, Ns, D_Li, t_ann, saturation_depth, r)

    push!(Ni_all, Ni_arr)
    push!(labels, "r = $(round(r*1e4, digits=1)) µm ")
end
# =========================
# PLOT Li distribution
# =========================
x_slices = collect(0:dx:Lx)
x_slices = x_slices[1:length(Ni_arr)]

x_slices = collect(0:dx:Lx)

p1 = plot(
    xlabel="depth (mm)",
    ylabel="N Li / slice",
    frame=:box,
    gridalpha=0.05,
    yticks=(0:5:maximum(Ni_all[1])),
    size=(400, 600),
    dpi=300
)

for i in 1:length(r_list)

    Ni_arr = Ni_all[i]
    x = x_slices[1:length(Ni_arr)] .* 10

    plot!(
        p1,
        x,
        Ni_arr,
        #yerr=sqrt.(Ni_arr),
        label=labels[i],
        lw=1.5,
    )
end

display(p1)
savefig(p1, "plot/max_li.png")
#=
# =========================
# CCE SIMULATION
# =========================
x_pos = 0:dx:0.11
N_charges = 10
N_repeat = 2

CCE_matrix = zeros(length(x_pos), N_repeat)

pbar = Progress(length(x_pos) * N_repeat)

for (i, x0) in enumerate(x_pos)
    for j in 1:N_repeat

        collected, trapped =
            multiple_charges_trapping_3D(
                x0, N_charges,
                cells, cell_size,
                saturation_depth, dx,
                Lx, Ly, Lz,
                Nt, σ,
                FCCD_cm
            )

        CCE_matrix[i, j] = collected / N_charges

        next!(pbar)
    end
end

CCE_mean = vec(mean(CCE_matrix, dims=2));
CCE_std = vec(std(CCE_matrix, dims=2));

# =========================
# PLOT CCE
# =========================
p2 = plot(
    x_pos .* 10,
    CCE_mean,
    ribbon=CCE_std,
    lw=2,
    xlabel="depth (mm)",
    ylabel="CCE",
    label="CCE",
    frame=:box
)

display(p2)

# =========================
# SAVE JSON
# =========================
filename = "JSON/Nd(x)-Nsat.json"

results = Dict(
    "x_pos" => collect(x_pos),
    "CCE_mean" => CCE_mean,
    "CCE_std" => CCE_std
)

if isfile(filename)
    println("File exists. Overwrite? (y/n)")
    if lowercase(readline()) != "y"
        println("Cancelled")
        return
    end
end

open(filename, "w") do f
    JSON.print(f, results, 4)
end

println("Saved: $filename")=#