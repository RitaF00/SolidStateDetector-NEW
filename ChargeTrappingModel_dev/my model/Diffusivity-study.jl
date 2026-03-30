
using Pkg
Pkg.activate(expanduser("~/Phd/SSDdev"))
Pkg.instantiate()
using Random
using Plots
using SpecialFunctions
using ProgressMeter
using Distributions
using JSON
using Unitful
using SolidStateDetectors
using SolidStateDetectors: Electron, Hole
T = Float64

#---------------
# test volume dimension
#---------------
Lx = 0.2
Ly = 0.2
Lz = 0.2

#---------------
# Li generation
#---------------
dx = 0.001   # Li diffusion slice in cm (10 µm)
α = 1.6 * 1e-11     # thinning factor
FCCD_cm = 0.1



#---------------
# Li-layer characteristics:
# - annealing time and tempertaures
# - donor concentration cm^-3
# - acceptor
# - neutral
#---------------
t_ann = 18.0 * 60.0
T_ann = 623.0
R = 1.98
H = (473 <= T_ann <= 873) ? 11800 : 10700
D0 = (473 <= T_ann <= 873) ? 2.5e-3 : 1.3e-3
Ns = 10^(21.27 - 2610 / T_ann)
D_Li = D0 * exp(-H / (R * T_ann))

Nd(x) = Ns .* erfc.(x ./ (2 * sqrt(D_Li * t_ann)))
N_neutral = 5.677e15  #cm-3
Na = 3e9 #cm^-3

#---------------
# Charge diffuusion
#---------------
Δt = 1.0      #ns
t_max = 10000 #ns
T_diff = 90.0 #k

# the function uses parameters for concentration expressed in m^-3
function _calculate_mobility_with_impurities_hole(
    Nn, bulk_imp, surface_imp, temperature)

    Ni = surface_imp .+ bulk_imp

    μI::T = 2.35e19 * temperature^1.5 / Ni / log(9.13e19 * temperature^2 / Ni) + 1.51e20 * temperature^1.5 / Ni / log(5.82e20 * temperature^2 / Ni)
    μA::T = 7.77e3 * temperature^-1.5
    μN::T = 1e2 / Nn * (2.31e18 + 2.36e20) * 0.82 * (0.228 * temperature^0.5 + 0.976 * temperature^-0.5)

    1 / (1 / μI + 1 / μA + 1 / μN)
end


# -----------------------------
# Generation cells with Li point-like
# -----------------------------
function generate_Li_cells(Lx, Ly, Lz, dx, α)

    cell_size = 0.0020 # dimension of the grid

    nx = Int(Lx / cell_size)
    ny = Int(Ly / cell_size)
    nz = Int(Lz / cell_size)

    cells = [Vector{NTuple{3,Float64}}() for _ in 1:nx, _ in 1:ny, _ in 1:nz]

    x_slices = 0:dx:Lx
    total_Li = 0

    for xi in x_slices

        Nd_val = Ns * erfc(xi / (2 * sqrt(D_Li * t_ann)))
        dV = dx * Ly * Lz
        λ = α * Nd_val * dV

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

    println("total number of Li generated = ", total_Li)

    return cells, nx, ny, nz, cell_size
end



# -----------------------------
# Trapping 
# -----------------------------
function multiple_charges_trapping_3D(x_charges, N, cells, cell_size)

    if isa(x_charges, Number)
        x_charges = fill(x_charges, N)
    end
    nx, ny, nz = size(cells)
    collected_count = 0

    for n in 1:N

        x = x_charges[n]
        y = rand() * Ly
        z = rand() * Lz
        Nd_data = Nd(x)
        # expressed in cm^2 V^-1 s^-1
        μh = _calculate_mobility_with_impurities_hole.(
            N_neutral .* 1e6,
            Na .* 1e6,
            Nd_data .* 1e6,
            T_diff) .* 1e4


        kB = 1.380649e-23   # J/K
        e = 1.602176634e-19 # C
        D = (kB * T_diff / e) * μh * 1e-9 # Conversion in cm²/ns
        σ = sqrt.(2 .* D .* Δt)

        for i in 1:Nt

            x += σ * randn()
            y += σ * randn()
            z += σ * randn()

            x = clamp(x, 0, Lx)
            y = clamp(y, 0, Ly)
            z = clamp(z, 0, Lz)

            ix = clamp(Int(floor(x / cell_size)) + 1, 1, nx)
            iy = clamp(Int(floor(y / cell_size)) + 1, 1, ny)
            iz = clamp(Int(floor(z / cell_size)) + 1, 1, nz)

            trapped = false

            for dxi in -1:1, dyi in -1:1, dzi in -1:1

                jx = ix + dxi
                jy = iy + dyi
                jz = iz + dzi

                if 1 ≤ jx ≤ nx && 1 ≤ jy ≤ ny && 1 ≤ jz ≤ nz

                    for (px, py, pz) in cells[jx, jy, jz]

                        if (x - px)^2 + (y - py)^2 + (z - pz)^2 < r_Li^2
                            trapped = true
                            break
                        end
                    end

                    if trapped
                        break
                    end
                end
            end

            if trapped
                break
            end

            if x == 0
                break
            end

            if x >= FCCD_cm
                collected_count += 1
                break
            end
        end
    end

    return collected_count
end


#---------- SIMUALTION --------------



# -----------------------------
# Li generation
# -----------------------------
cells, nx, ny, nz, cell_size = generate_Li_cells(Lx, Ly, Lz, dx, α)
n_precipitates = sum(length(cells[ix, iy, iz]) for ix in 1:nx, iy in 1:ny, iz in 1:nz)
println("Numero precipitati Li = ", n_precipitates)

# -----------------------------
# CCE curve
# -----------------------------
x_pos = 0:0.002:0.11

N_charges = 500
N_repeat = 25

N_matrix = zeros(Int, length(x_pos), N_repeat)

pbar = Progress(length(x_pos) * N_repeat)

for (i, x0) in enumerate(x_pos)
    for j in 1:N_repeat

        N_matrix[i, j] = multiple_charges_trapping_3D(x0, N_charges, cells, cell_size)

        next!(pbar)
    end
end


CCE = N_matrix ./ N_charges
CCE_mean = vec(mean(CCE, dims=2))
CCE_std = vec(std(CCE, dims=2))

# -----------------------------
# Plot CCE
# -----------------------------
plt = plot(x_pos, CCE_mean,
    lw=2,
    color=:orange,
    xlabel="depth (cm)",
    ylabel="CCE",
    label="20 μm precipitates",
    size=(800, 600),
    frame=:box
)

plot!(
    x_pos, CCE_mean,
    yerr=CCE_std,
    seriestype=:scatter,
    ms=2,
    label=false
)

vline!(plt, [FCCD_cm], linestyle=:dash, label="FCCD")

savefig(plt, "CCE_diffusivity.png")
display(plt)




#----dizionario---
results = Dict(
    "x_pos" => collect(x_pos),        # range → array
    "CCE_mean" => CCE_mean,
    "CCE_std" => CCE_std
)

# -----------------------------
# Salvataggio su file JSON
# -----------------------------
open("Diffusivity_CCE.json", "w") do f
    JSON.print(f, results, 4)  # 4 = indentazione bella leggibile
end

println("File JSON salvato: Diffusivity_CCE.json")






