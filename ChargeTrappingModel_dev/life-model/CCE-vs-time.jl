using Random
using Plots
using SpecialFunctions
using Statistics
using Base.Threads
using ProgressMeter

gr()

# ============================================================
# GEOMETRY (cm)
# ============================================================
const Lx = 0.2
const Ly = 0.2
const Lz = 0.2

const dx = 0.001
#const FCCD_cm = 0.1
const FCCD_cm = 0.08957282
const x_RDR = 0.065

# ============================================================
# PHYSICS
# ============================================================
const t_ann = 18 * 60
const T_ann = 623

const R = 1.98
const H = 11800
const D0 = 2.5e-3

const Ns = 10^(21.27 - 2610 / T_ann)
const D_Li = D0 * exp(-H / (R * T_ann))

# hole diffusion
const D = 28.9
const Δt = 1.0          # ns
const Nt_steps = 10000

const T_diff = 90

# Li trap radius
const r_Li = 0.002

# ============================================================
# PRECOMPUTED CONSTANTS
# ============================================================
const m0 = 9.11e-31
const m_eff = 0.21 * m0
const kB = 1.38e-23

const v_th = sqrt(3 * kB * T_diff / m_eff) * 100
const σ_trap = π * r_Li^2

const pref_τ = 1.0 / (v_th * σ_trap)

# diffusion step (cm)
const step_sigma = sqrt(2 * D * Δt) * 1e-4

# ============================================================
# LUT
# ============================================================
const LUT_dx = 1e-4
const xgrid = 0:LUT_dx:FCCD_cm

# ============================================================
# SIMULATION
# ============================================================
step = 0.005
x_pos = step:step:FCCD_cm-step

x_mm = x_pos .* 10

N_charges = 100
N_repeat = 50

Nx = length(x_pos)

# ============================================================
# α VALUES
# ============================================================
alphas = [
    0,
    #1e-12,
    #1e-11,
    #1e-10,
    #1e-9,
    1e-8,
]

# ============================================================
# RNG
# ============================================================
rngs = [
    MersenneTwister(1234 + i)
    for i in 1:nthreads()
]

# ============================================================
# τ(x)
# ============================================================
function τ_exact(x, α)

    Nd = Ns * erfc(
        x / (2 * sqrt(D_Li * t_ann))
    )

    return pref_τ / (α * Nd)
end

function build_tau_grid(alpha)

    return [
        τ_exact(x, alpha)
        for x in xgrid
    ]
end

@inline function τ_hole(x, τgrid)

    idx = Int(clamp(
        floor(x / LUT_dx) + 1,
        1,
        length(τgrid)
    ))

    return τgrid[idx]
end

# ============================================================
# trapping probability
# ============================================================
@inline function trapping_probability(x, τgrid)

    # no trapping outside Li region
    if x >= x_RDR
        return 0.0
    end

    τ_ns = τ_hole(x, τgrid) * 1e9

    if τ_ns <= 0
        return 1.0
    end

    return 1 - exp(-Δt / τ_ns)
end

# ============================================================
# propagate charge
# ============================================================
@inline function propagate_charge(x0, rng, τgrid)

    x = x0
    y = Ly / 2
    z = Lz / 2

    # surface loss immediata
    if x <= dx
        return (:surface_loss, 0.0)
    end

    for i in 1:Nt_steps

        x += step_sigma * randn(rng)
        y += step_sigma * randn(rng)
        z += step_sigma * randn(rng)

        y = clamp(y, 0, Ly)
        z = clamp(z, 0, Lz)

        # surface loss
        if x <= 0.0
            return (:surface_loss, i * Δt)
        end

        # trapping
        if rand(rng) < trapping_probability(x, τgrid)
            return (:trapped, i * Δt)
        end

        # collected
        if x >= FCCD_cm
            return (:collected, i * Δt)
        end
    end

    return (:trapped, Nt_steps * Δt)
end

# ============================================================
# run depth
# ============================================================
function run_depth(x0, N, rng, τgrid)

    collected = 0
    times = Float64[]

    for i in 1:N

        status, t = propagate_charge(x0, rng, τgrid)

        push!(times, t)

        if status == :collected
            collected += 1
        end
    end

    return collected, times
end

# ============================================================
# PLOTS
# ============================================================
p1 = plot(
    xlabel="Depth (mm)",
    ylabel="CCE",
    frame=:box,
    lw=2,
    legend=:bottomright
)

p2 = plot(
    xlabel="Depth (mm)",
    ylabel="t₉₀ - t₁₀ (μs)",
    frame=:box,
    lw=2,
    legend=:topright
)

p3 = plot(
    xlabel="CCE",
    ylabel="t₉₀ - t₁₀ (μs)",
    frame=:box,
    lw=2,
    legend=:topright
)

p4 = plot(
    xlabel="Depth (mm)",
    ylabel="t₉₀₋₁₀(α) / t₉₀₋₁₀(0)",
    frame=:box,
    lw=2,
    legend=:topright
)

# ============================================================
# STORAGE
# ============================================================
reference_timing = Vector{Float64}()

heatmaps = []

# ============================================================
# GLOBAL PROGRESS
# ============================================================
total_steps = length(alphas) * Nx * N_repeat

pbar = Progress(total_steps)

# ============================================================
# MAIN LOOP
# ============================================================
for (ia, alpha_val) in enumerate(alphas)

    println("\nRunning α = $alpha_val")

    τgrid = build_tau_grid(alpha_val)

    times_by_depth =
        Vector{Vector{Vector{Float64}}}(undef, Nx)

    N_matrix = zeros(Int, Nx, N_repeat)

    # ========================================================
    # PARALLEL
    # ========================================================
    Threads.@threads for i in eachindex(x_pos)

        rng = rngs[threadid()]

        x0 = x_pos[i]

        rep_times =
            Vector{Vector{Float64}}(undef, N_repeat)

        for j in 1:N_repeat

            collected, times =
                run_depth(x0, N_charges, rng, τgrid)

            N_matrix[i, j] = collected

            rep_times[j] = times

            next!(pbar)
        end

        times_by_depth[i] = rep_times
    end

    # ========================================================
    # CCE
    # ========================================================
    CCE = N_matrix ./ N_charges

    CCE_mean = vec(mean(CCE, dims=2))
    CCE_std = vec(std(CCE, dims=2))

    plot!(
        p1,
        x_mm,
        CCE_mean,
        ribbon=CCE_std,
        fillalpha=0.25,
        label="α = $(alpha_val)"
    )

    # ========================================================
    # t90-t10
    # ========================================================
    mean_dt = Float64[]
    std_dt = Float64[]

    for i in eachindex(x_pos)

        reps = times_by_depth[i]

        dt_rep = Float64[]

        for rep in reps

            if length(rep) > 10

                t10 = quantile(rep, 0.10)
                t90 = quantile(rep, 0.90)

                push!(
                    dt_rep,
                    (t90 - t10) * 0.001
                )
            end
        end

        push!(mean_dt, mean(dt_rep))
        push!(std_dt, std(dt_rep))
    end

    plot!(
        p2,
        x_mm,
        mean_dt,
        ribbon=std_dt,
        fillalpha=0.25,
        label="α = $(alpha_val)"
    )

    # ========================================================
    # t(CCE)
    # ========================================================
    plot!(
        p3,
        CCE_mean,
        mean_dt,
        lw=2,
        marker=:circle,
        label="α = $(alpha_val)"
    )

    # ========================================================
    # normalized timing
    # ========================================================
    if ia == 1

        global reference_timing
        reference_timing = copy(mean_dt)

    else

        ratio = mean_dt ./ reference_timing

        plot!(
            p4,
            x_mm,
            ratio,
            lw=2,
            marker=:circle,
            label="α = $(alpha_val)"
        )
    end

    # ========================================================
    # HEATMAP
    # ========================================================
    all_x = Float64[]
    all_t = Float64[]

    for (i, x0) in enumerate(x_pos)

        reps = times_by_depth[i]

        for rep in reps

            append!(
                all_x,
                fill(x0 * 10, length(rep))
            )

            append!(
                all_t,
                rep .* 0.001
            )
        end
    end

    # ============================================================
    # HEATMAP MIGLIORATA
    # ============================================================

    dx_mm = step * 10
    x_bins = 0:0.075:0.9

    dt_us = Δt * 0.01
    t_bins = 0:dt_us:maximum(all_t)

    h = histogram2d(
        all_x,
        all_t,
        bins=(x_bins, t_bins),
        xlabel="Depth (mm)",
        ylabel="Collection time (μs)",
        colorbar_title="counts",
        frame=:box
    )

    h = histogram2d(
        all_x,
        all_t,
        bins=(x_bins, t_bins),
        xlabel="Depth (mm)",
        ylabel="Collection time (μs)",
        color=:plasma,
        title="α = $(alpha_val)",
        colorbar_title="counts",
        frame=:box,
        xlim=(0, 0.9),
        ylim=(minimum(all_t), maximum(all_t)),
        normalize=:none
    )

    push!(heatmaps, h)
end

# ============================================================
# FINAL FIGURES
# ============================================================
pfinal1 = plot(
    p1,
    p2,
    p3,
    p4,
    layout=(2, 2),
    size=(1100, 900),
    dpi=300
)

savefig(
    pfinal1,
    "CCE_timing_diagnostics.png"
)

# ============================================================
# HEATMAP PANEL
# ============================================================
nH = length(heatmaps)

rows = ceil(Int, nH / 2)

pfinal2 = plot(
    heatmaps...,
    layout=(rows, 2),
    size=(1200, 300 * rows),
    dpi=300
)

savefig(
    pfinal2,
    "timing_heatmaps.png"
)

println("\nDONE")
println("Saved:")
println("CCE_timing_diagnostics.png")
println("timing_heatmaps.png")