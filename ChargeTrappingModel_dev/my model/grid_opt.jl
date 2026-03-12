using Random
using Plots
using SpecialFunctions
using ProgressMeter
using Distributions

# -----------------------------
# Parametri fisici e griglia
# -----------------------------
Lx = 0.2      # profondità in cm
Ly = 0.2       # estensione y in cm
Lz = 0.2       # estensione y in cm

dx = 0.0001     # passo x per generare slice 1 μm
α = 1e-9       # thinning factor per Poisson

FCCD_cm = 0.09
PN_cm = 0.11

#annealing
t_ann = 18 * 60       # 18 min in sec
T_ann = 623

#diffusione
R = 1.98
H = 11800
D0 = 2.5e-3
Ns = 10^(21.27 - 2610 / T_ann)
D_Li = D0 * exp(-H / (R * T_ann))

# Parametri diffusione
D = 28.9           # diffusione [µm^2/ns]
Δt = 1.0           # ns
t_max = 10000      # ns
Nt = Int(t_max / Δt)
σ = sqrt(2 * D * Δt) * 1e-4  # µm -> cm

# Raggio per cattura Li (ora corrisponde a un voxel)
r_Li = 0.002  # cm (20 µm)

# -----------------------------
# Funzione per generare matrice di Li
# -----------------------------
function generate_Li_grid(Lx, Ly, Lz, dx, α)
    r_Li = 0.002
    nx = Int(Lx / r_Li)
    ny = Int(Ly / r_Li)
    nz = Int(Lz / r_Li)
    Li_grid = falses(nx, ny, nz)

    x_slices = 0:dx:Lx
    for xi in x_slices
        Nd_val = Ns * erfc(xi / (2 * sqrt(D_Li * t_ann)))
        #c0 = 0.005
        #L_max = 0.06
        #Nd_val = c0 * max(0, 1 - xi / L_max)
        dV = dx * Ly * Lz
        λ = α * Nd_val * dV
        N_i = rand(Poisson(λ))
        for _ in 1:N_i
            x_pos = xi + rand() * dx
            y_pos = rand() * Ly
            z_pos = rand() * Lz
            ix = clamp(Int(floor(x_pos / r_Li)) + 1, 1, nx)
            iy = clamp(Int(floor(y_pos / r_Li)) + 1, 1, ny)
            iz = clamp(Int(floor(z_pos / r_Li)) + 1, 1, nz)
            Li_grid[ix, iy, iz] = true
        end
    end

    return Li_grid, nx, ny, nz
end

function multiple_charges_trapping_3D(x_charges, N::Int=100, Li_grid=nothing, Lx=Lx, Ly=Ly, Lz=Lz)
    if isa(x_charges, Number)
        x_charges = fill(x_charges, N)
    elseif length(x_charges) != N
        error("lunghezza x_charges != N")
    end

    nx, ny, nz = size(Li_grid)
    r_Li = Lx / nx

    collected_count = 0
    for n in 1:N
        x_charge = x_charges[n]
        y_charge = rand() * Ly
        z_charge = rand() * Lz
        trapped = false

        if x_charge >= Lx
            collected_count += 1
            continue
        end

        for i in 1:Nt
            x_charge += σ * randn()
            y_charge += σ * randn()
            z_charge += σ * randn()

            x_charge = clamp(x_charge, 0, Lx)
            y_charge = clamp(y_charge, 0, Ly)
            z_charge = clamp(z_charge, 0, Lz)

            ix = clamp(Int(floor(x_charge / r_Li)) + 1, 1, nx)
            iy = clamp(Int(floor(y_charge / r_Li)) + 1, 1, ny)
            iz = clamp(Int(floor(z_charge / r_Li)) + 1, 1, nz)

            if Li_grid[ix, iy, iz]
                trapped = true
                break
            end

            if x_charge >= FCCD_cm
                collected_count += 1
                break
            end
        end
    end

    return collected_count
end

# -----------------------------
# Generazione matrice Li 3D
# -----------------------------
Li_grid, nx, ny, nz = generate_Li_grid(Lx, Ly, Lz, dx, α)
println("Matrice Li 3D generata!")


factor = 5  # riduce il numero di voxel da visualizzare
Li_plot = Li_grid[1:factor:end, 1:factor:end, 1:factor:end]

nx, ny, nz = size(Li_plot)
xs, ys, zs = Float64[], Float64[], Float64[]

for ix in 1:nx
    for iy in 1:ny
        for iz in 1:nz
            if Li_plot[ix, iy, iz]
                push!(xs, (ix - 0.5) * dx * factor)  # coord in cm
                push!(ys, (iy - 0.5) * dx * factor)
                push!(zs, (iz - 0.5) * dx * factor)
            end
        end
    end
end

p = scatter3d(xs, ys, zs, markersize=3, markercolor=:red,
    xlabel="x (cm)", ylabel="y (cm)", zlabel="z (cm)",
    title="Distribuzione 3D dei precipitatii Li",
    xlim=(0, 0.2),
    ylim=(0, 0.2),
    zlim=(0, 0.2))
#display(p)
# -----------------------------
# Simulazione CCE in funzione della profondità
# -----------------------------
x_pos = 0:0.001:0.15
N_charges = 10000
N_repeat = 5

N_matrix = zeros(Int, length(x_pos), N_repeat)
p = Progress(length(x_pos) * N_repeat, desc="Simulazione CCE")

for (i, x0) in enumerate(x_pos)
    for j in 1:N_repeat
        N_matrix[i, j] = multiple_charges_trapping_3D(x0, N_charges, Li_grid)
        next!(p)
    end
end

CCE = N_matrix ./ N_charges
CCE_mean = vec(mean(CCE, dims=2))
CCE_std = vec(std(CCE, dims=2))

# -----------------------------
# Plot
# -----------------------------
xticks_positions = 0:0.02:0.15

plt = plot(x_pos, CCE_mean, ribbon=CCE_std, color=:orange,
    xlabel="depth (cm)", ylabel="CCE", legend=false, xticks=(xticks_positions, string.(round.(xticks_positions, digits=2))))
vline!(plt, [PN_cm], ls=:dash, color=:green, label="P-N junction")
display(plt)

