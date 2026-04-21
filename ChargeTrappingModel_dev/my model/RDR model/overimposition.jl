using Random
using Plots
using SpecialFunctions
using Distributions
using ProgressMeter
using Printf

gr()

"""
In questo codice vedo quante litio posso mettere prima di saturare la fetta/sistema
restituisce solo la percentuale dei litio posizionati.

se la posizione del litio è già occupata da un precipitato fissato in precedenza,
il sistema prova fino a 250 volte a piazzarlo da un'altra parte.
Se non riesce, il litio viene buttato.

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

function run_simulation(r, Nrun)

    # -------------------------
    # DOMAIN
    # -------------------------
    Lx = 0.2
    Ly = 0.2
    Lz = 0.2
    dx = 0.001

    α = 1.6e-11
    Nd_saturation = 1e14

    # -------------------------
    # PHYSICS
    # -------------------------
    t_ann = 18 * 60
    T_ann = 623

    D0 = 2.5e-3
    H = 11800
    R = 1.98

    Ns = 10^(21.27 - 2610 / T_ann)
    D_Li = D0 * exp(-H / (R * T_ann))

    depth_list = 0:dx:0.11
    Nd_vals = Ns .* erfc.(depth_list ./ (2 * sqrt(D_Li * t_ann)))
    idx = findfirst(<=(0), Nd_vals .- Nd_saturation)
    saturation_depth = depth_list[idx]

    # -------------------------
    # STORAGE
    # -------------------------
    results = Float64[]
    failures = Int[]

    V_particle = (4 / 3) * π * r^3
    V_slice = dx * Ly * Lz

    # -------------------------
    # RUNS
    # -------------------------
    for run in 1:Nrun

        failed = 0
        total_accepted = 0

        cell_size = r
        nx = Int(floor(Lx / cell_size))
        ny = Int(floor(Ly / cell_size))
        nz = Int(floor(Lz / cell_size))

        cells = [Vector{LiParticle}() for _ in 1:nx, _ in 1:ny, _ in 1:nz]

        Ni_arr = zeros(Int, Int(floor(Lx / dx)))

        for (i, xi) in enumerate(0:dx:Lx)

            if xi >= saturation_depth
                continue
            end

            Nd_val = Ns * erfc(xi / (2 * sqrt(D_Li * t_ann)))
            λ = α * max(Nd_val - Nd_saturation, 0.0) * dx * Ly * Lz

            N_i = rand(Poisson(λ))

            for _ in 1:N_i

                placed = false

                for _ in 1:250

                    x_pos = xi + rand() * dx
                    y_pos = rand() * Ly
                    z_pos = rand() * Lz

                    ix = clamp(Int(floor(x_pos / cell_size)) + 1, 1, nx)
                    iy = clamp(Int(floor(y_pos / cell_size)) + 1, 1, ny)
                    iz = clamp(Int(floor(z_pos / cell_size)) + 1, 1, nz)

                    overlap = false

                    for dxs in -1:1, dys in -1:1, dzs in -1:1
                        jx, jy, jz = ix + dxs, iy + dys, iz + dzs

                        if 1 ≤ jx ≤ nx && 1 ≤ jy ≤ ny && 1 ≤ jz ≤ nz
                            for p in cells[jx, jy, jz]
                                if (x_pos - p.x)^2 + (y_pos - p.y)^2 + (z_pos - p.z)^2 < (2r)^2
                                    overlap = true
                                    break
                                end
                            end
                        end

                        if overlap
                            break
                        end
                    end

                    if !overlap
                        push!(cells[ix, iy, iz], LiParticle(x_pos, y_pos, z_pos, r))
                        Ni_arr[i] += 1
                        placed = true
                        total_accepted += 1
                        break
                    end
                end

                if !placed
                    failed += 1
                end
            end
        end

        push!(results, sum(Ni_arr))
        push!(failures, failed)
    end

    # -------------------------
    # GLOBAL METRICS
    # -------------------------
    total_accepted = sum(results)
    total_failed = sum(failures)

    survival = total_accepted + total_failed > 0 ?
               total_accepted / (total_accepted + total_failed) : 0.0

    mean_particles = mean(results)
    std_particles = std(results)

    packing_mean = mean(results) * V_particle / (Lx * Ly * Lz)
    packing_std = std(results) * V_particle / (Lx * Ly * Lz)

    # -------------------------
    # FINAL OUTPUT
    # -------------------------


    println("\n===== FINAL BEHAVIOR =====")
    @printf("radius r = %.2g µm\n", r * 1e4)
    println("runs = $Nrun")
    #=
    println("\n--- PARTICLE STATISTICS ---")
    @printf("mean accepted particles = %.2g\n", mean_particles)
    @printf("std accepted particles = %.2g\n", std_particles)

    println("\n--- SYSTEM EFFICIENCY ---")
    =#
    println("total produced = $(total_failed+total_accepted)")
    println("total accepted = $total_accepted")
    println("total failed = $total_failed")
    @printf("global survival rate = %.2g %%\n", survival * 100)

    println("\n--- PACKING ---")
    @printf("mean packing fraction = %.2g\n", packing_mean)
    @printf("packing std = %.2g\n", packing_std)

    return results, failures, survival
end

# =========================
# RUN
# =========================
r = 0.012   # 20 µm
Nrun = 25

results, failures, sr = run_simulation(r, Nrun)