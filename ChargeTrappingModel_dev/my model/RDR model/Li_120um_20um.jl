using Random
using Plots
using SpecialFunctions
using Distributions
using JSON
using ProgressMeter
gr()

"""
In questo codice utilizzo due differenti valori di raggio
proprio per questo modtivo devo tenere in considerazione quanto ne posso avere

Il numero di litio prodotto per cella viene  preso come il minimo tra il numero di eventi geoemtrico (V_slice/V_litio) 
e il numero di litio ottentuto tramite la poissoniana.

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
# PARAMETERS
# =========================
Lx = 0.2
Ly = 0.2
Lz = 0.2

dx = 0.001

α = 1.6e-11
Nd_saturation = 1e14

r_Li = 0.002

# =========================
# ANNEALING
# =========================
t_ann = 18 * 60
T_ann = 623
D0 = 2.5e-3
H = 11800
R = 1.98

Ns = 10^(21.27 - 2610 / T_ann)
D_Li = D0 * exp(-H / (R * T_ann))

# =========================
# SATURATION DEPTH
# =========================
depth_list = 0:dx:0.11
Nd_vals = Ns .* erfc.(depth_list ./ (2 * sqrt(D_Li * t_ann)))

idx = findfirst(i -> Nd_vals[i] - Nd_saturation <= 0, eachindex(depth_list))
saturation_depth = depth_list[idx]

println("Saturation depth = ", saturation_depth * 10, " mm")

# =========================
# RADIUS PROFILE
# =========================
function r_Li_profile(x)
    return x <= saturation_depth / 2 ? 6 * r_Li : r_Li
end
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
# OVERLAP CHECK (FAST CELL LIST)
# =========================
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

# =========================
# GENERATION (PHYSICALLY CONSISTENT)
# =========================
function generate_Li()

    cell_size = 2 * r_Li
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

        # -------------------------
        # local physical density
        # -------------------------
        Nd_val = Ns * erfc(xi / (2 * sqrt(D_Li * t_ann)))
        density = α * max(Nd_val - Nd_saturation, 0.0)

        λ = density * dx * Ly * Lz
        r = r_Li_profile(xi)
        # geometric upper bound (VERY IMPORTANT)
        V_particle = (4 / 3) * π * r^3
        V_slice = dx * Ly * Lz
        N_max = Int(floor(V_slice / V_particle))
        N_target = min(Int(round(λ)), N_max)
        N_target = truncated_poisson(λ, N_max)

        Ni_geom[i] = N_target

        accepted = 0

        # -------------------------
        # deposition WITHOUT retry
        # -------------------------
        trials = N_target * 3

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

        #println("slice $i | target=$N_target | accepted=$accepted")
    end

    println("\nTotal accepted particles = $total")

    return cells, Ni_geom, Ni_accept, x_slices
end

# =========================
# RUN
# =========================
cells, Ni_geom, Ni_accept, x_slices = generate_Li()

# =========================
# EXTRACT POINTS
# =========================
xs, ys, zs, sizes, colors = Float64[], Float64[], Float64[], Float64[], Any[]

for cell in cells
    for p in cell
        push!(xs, p.x)
        push!(ys, p.y)
        push!(zs, p.z)

        if p.x <= saturation_depth / 2
            push!(sizes, 6)
            push!(colors, :blue)
        else
            push!(sizes, 2)
            push!(colors, :red)
        end
    end
end

# =========================
# 3D PLOT
# =========================
p3d = scatter3d(xs, ys, zs,
    markersize=sizes,
    markercolor=colors,
    alpha=0.7,
    xlabel="x",
    ylabel="y",
    zlabel="z",
    legend=false,
    aspect_ratio=:equal,
    size=(800, 800)
)

display(p3d)

# =========================
# PROJECTIONS
# =========================
p_xy = scatter(xs, ys,
    markersize=sizes,
    markercolor=colors,
    xlabel="x",
    ylabel="y",
    legend=false,
    title="XY",
    frame=:box
)


p_xz = scatter(xs, zs,
    markersize=sizes,
    markercolor=colors,
    xlabel="x",
    ylabel="z",
    legend=false,
    title="XZ",
    frame=:box
)

p_yz = scatter(ys, zs,
    markersize=sizes,
    markercolor=colors,
    xlabel="y",
    ylabel="z",
    legend=false,
    title="YZ",
    frame=:box
)

display(plot(p_xy, p_xz, p_yz, layout=(1, 3), size=(2000, 600)))


# =========================
# SLICE PLOT (IMPORTANT FIX)
# =========================
p = plot()
vspan!(p, [0, saturation_depth * 10 / 2],
    color=:deepskyblue,
    alpha=0.2,
    label=""
)
plot!(p, x_slices .* 10, Ni_geom,
    lw=2, color=:orange, label="",
    xlim=[0, 0.7])
xlabel!("x (mm)")
ylabel!("Li per slice")
plot!(legend=:topright, frame=:box, size=(700, 500))
vline!([saturation_depth * 10 / 2], ls=:dash, color=:grey, alpha=0.25, label="RDR/2")

display(p)