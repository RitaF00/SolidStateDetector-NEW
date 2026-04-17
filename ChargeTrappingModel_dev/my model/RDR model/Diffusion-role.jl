using Random
using Plots
using SpecialFunctions
using Distributions
using ProgressMeter
using Printf
gr()

"""
In questo codice vedo il ruolo della diffusione sulla curva CCE
ATTENZIONE : ricoda che quando cambi la dimensione del precipitato di litio devi anche
cambaire la dimensione della cella altrimenti SOVRASTIMI LA CCE
"""

# =====================
# Li generation (FIXED)
# =====================
function generate_Li_cells(Lx, Ly, Lz, dx, α, Ns, D_Li, t_ann, x_saturation, Nd_saturation, r_Li)

    cell_size = r_Li  # cm
    cell_size = 0.02

    nx = Int(floor(Lx / cell_size))
    ny = Int(floor(Ly / cell_size))
    nz = Int(floor(Lz / cell_size))

    cells = [Vector{NTuple{3,Float64}}() for _ in 1:nx, _ in 1:ny, _ in 1:nz]

    x_slices = 0:dx:Lx
    total_Li = 0

    for xi in x_slices

        if xi >= x_saturation
            continue
        end

        Nd_val = Ns * erfc(xi / (2 * sqrt(D_Li * t_ann)))

        # CORREZIONE FISICA
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
    return cells, nx, ny, nz, cell_size, total_Li
end


# =====================
# Charge transport
# =====================
function multiple_charges_trapping_3D(
    x_charges, N, cells, cell_size,
    x_saturation, dx, Lx, Ly, Lz,
    Nt, σ, FCCD_cm, r_Li
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

            x += σ * randn()
            y += σ * randn()
            z += σ * randn()

            x = clamp(x, 0, Lx)
            y = clamp(y, 0, Ly)
            z = clamp(z, 0, Lz)

            if x >= FCCD_cm
                collected_count += 1
                break
            end

            if x <= dx
                trapped = true
                break
            end

            if x <= x_saturation

                ix = clamp(Int(floor(x / cell_size)) + 1, 1, nx)
                iy = clamp(Int(floor(y / cell_size)) + 1, 1, ny)
                iz = clamp(Int(floor(z / cell_size)) + 1, 1, nz)

                for dxi in -1:1, dyi in -1:1, dzi in -1:1

                    jx, jy, jz = ix + dxi, iy + dyi, iz + dzi

                    if 1 ≤ jx ≤ nx && 1 ≤ jy ≤ ny && 1 ≤ jz ≤ nz

                        for (px, py, pz) in cells[jx, jy, jz]

                            if (x - px)^2 + (y - py)^2 + (z - pz)^2 < r_Li^2
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
Lx, Ly, Lz = 0.2, 0.2, 0.2
dx = 0.001
FCCD_cm = 0.1

D = 28.9 # μm^2 / ns
Δt = 1.0
t_max = 10000
Nt = Int(t_max / Δt)

D_arr = [10 * D, 5 * D, 3 * D, D]


σ_arr = sqrt.(6 .* D_arr .* Δt) * 1e-4 #cm

println("============================")
println("Diffusion $(σ_arr.*1e4)")
println("============================")

r_Li = 0.002
r_Li = 0.006

t_ann = 18 * 60
T_ann = 623

R = 1.98
H = 11800
D0 = 2.5e-3

Ns = 10^(21.27 - 2610 / T_ann)
D_Li = D0 * exp(-H / (R * T_ann))

α = 1.6e-11


# =====================
# SCAN Nd_saturation
# =====================
Nd_sat_list = [1e17, 1e16, 1e15, 1e13]
colors = [:mediumpurple2, :deeppink, :deepskyblue, :orange, :lime]

x_pos = 0:dx:0.11
N_charges = 500
N_repeat = 25

CCE_all = Dict()
precip_counts = Dict()


# =====================
# MAIN LOOP
# =====================
for (k, σ) in enumerate(σ_arr)
    local depth_list, Nd_vals, Nd_saturation, diff_vals, idx,
    saturation_depth, cells, nx, ny, nz, cell_size,
    CCE_matrix, pbar, CCE_mean, CCE_std

    println("===============    σ = $( @sprintf("%.2g", σ*1e4) ) μm   ===============")

    depth_list = 0:dx:0.11
    Nd_vals = Ns .* erfc.(depth_list ./ (2 * sqrt(D_Li * t_ann)))

    Nd_saturation = 1e14

    diff_vals = Nd_vals .- Nd_saturation
    idx = findfirst(i -> diff_vals[i] <= 0, eachindex(diff_vals))

    if idx === nothing
        println("No saturation reached")
        continue
    end

    saturation_depth = depth_list[idx]
    print("saturation depth $saturation_depth")

    cells, nx, ny, nz, cell_size, total_Li =
        generate_Li_cells(
            Lx, Ly, Lz, dx, α, Ns, D_Li, t_ann,
            saturation_depth, Nd_saturation, r_Li
        )

    precip_counts[Nd_saturation] = total_Li


    # RAW STORAGE
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
                    FCCD_cm, r_Li
                )

            CCE_matrix[i, j] = collected / N_charges
            next!(pbar)
        end
    end


    # NORMALIZATION + STATS


    CCE_mean = vec(mean(CCE_matrix, dims=2))
    CCE_std = vec(std(CCE_matrix, dims=2))

    CCE_all[σ] = (CCE_mean, CCE_std)
end


# =====================
# PLOT
# =====================
p = plot(
    xlabel="depth (mm)",
    ylabel="CCE",
    frame=:box,
    title="Li size = $(r_Li*1e4) μm",
    dpi=300,
    size=(450, 600)
)

for (k, σ) in enumerate(σ_arr)

    if haskey(CCE_all, σ)

        local CCE_mean, CCE_std

        CCE_mean, CCE_std = CCE_all[σ]

        plot!(
            p,
            x_pos .* 10,
            CCE_mean,
            lw=2,
            color=colors[k],
            label="σ = $( @sprintf("%.2g", σ*1e4) ) μm",
            ribbon=CCE_std,
            fillalpha=0.25
        )
    end
end

display(p)
savefig(p, "plot/CCE_Diffusion_Rli=$(r_Li*1e4)_μm_bis2.png")


# =====================
# SAVE JSON
# =====================
filename = "JSON/CCE_Diffusion_Rli=$(r_Li*1e4)_μm_bis2.json"

results = Dict()

for (σ, (mean_vals, std_vals)) in CCE_all
    key = "sigma_$(round(σ*1e4, digits=2))um"

    results[key] = Dict(
        "CCE_mean" => mean_vals,
        "CCE_std" => std_vals
    )
end

results["x_pos"] = collect(x_pos)

if isfile(filename)
    println("Il file $filename esiste già. Vuoi sovrascriverlo? (y/n)")
    risposta = readline()

    if lowercase(risposta) != "y"
        println("Operazione annullata.")
    else
        open(filename, "w") do f
            JSON.print(f, results, 4)
        end
        println("File JSON salvato: $filename")
    end
else
    open(filename, "w") do f
        JSON.print(f, results, 4)
    end
    println("File JSON salvato: $filename")
end