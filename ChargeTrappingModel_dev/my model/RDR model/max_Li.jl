using Random
using Plots
using SpecialFunctions
using Distributions
using JSON
using ProgressMeter

gr()

"""
In questo codice guardo qual è il nuemro massimo di precipittai che posso mettere
per ogni fetta facendo variare il raggio.
"""

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

    for _ in 1:100
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

        f = 1

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

#=
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
    yticks=(0:1:maximum(Ni_all[1])),
    ylim=[0, 10],
    size=(400, 600),
    #dpi=300
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
=#


"""
produzione dei precipitati utilizzando come volume 100 volte il raggio
"""


function generate_Li_cells_prop_r_size(dx, α,
    Ns, D_Li, t_ann,
    x_saturation, r
)
    Lx = Ly = Lz = 100 * r
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

        f = 1

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




α = 1.6e-11
alpha_list = [1.6e-11, 2e-12, 5e-13]
alpha_list = [1.6e-11, 1.6e-11, 1.6e-11]
r_list = [0.002, 0.006, 0.012]

Ni_all = []
labels = []
for (i, r) in enumerate(r_list)
    println("r = $(r*1e4) μm")
    global cells, nx, ny, nz, cell_size, Ni_arr
    cells, nx, ny, nz, cell_size, Ni_arr =
        generate_Li_cells_prop_r_size(dx, alpha_list[i], Ns, D_Li, t_ann, saturation_depth, r)

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
    yticks=(0:50:maximum(Ni_all[2])),
    size=(400, 600),
    #dpi=300
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
savefig(p1, "plot/max_li_volumesize=100r.png")

println(sum(Ni_all[1]))
println(sum(Ni_all[2]))
println(sum(Ni_all[3]))