using Random
using Plots
using SpecialFunctions
using Statistics
using Base.Threads
using ProgressMeter

gr()

# ============================================================
# GEOMETRIA (cm)
# ============================================================
const Lx = 0.2
const Ly = 0.2
const Lz = 0.2

const dx = 0.001
const FCCD_cm = 0.1
const x_RDR = 0.065

# ============================================================
# PARAMETRI FISICI
# ============================================================
const t_ann = 18 * 60
const T_ann = 623

const R = 1.98
const H = 11800
const D0 = 2.5e-3

const Ns = 10^(21.27 - 2610 / T_ann)
const D_Li = D0 * exp(-H / (R * T_ann))

const D = 28.9
const Δt = 1.0
const Nt = 10000
const T_diff = 90

const r_Li = 0.002

# ============================================================
# COSTANTI PRECALCOLATE
# ============================================================
const m0 = 9.11e-31
const m_eff = 0.21 * m0
const kB = 1.38e-23

const v_th = sqrt(3 * kB * T_diff / m_eff) * 100
const σ_trap = π * r_Li^2

const pref_τ = 1.0 / (v_th * σ_trap)

const step_sigma = sqrt(2 * D * Δt) * 1e-4

# ============================================================
# LOOKUP TABLE
# ============================================================
const LUT_dx = 1e-4
const xgrid = 0:LUT_dx:FCCD_cm

# ============================================================
# SIMULATION PARAMETERS
# ============================================================
step = 0.005
x_pos = step:step:FCCD_cm-step

x_mm = x_pos .* 10

N_charges = 1000
N_repeat = 1000

Nx = length(x_pos)

xticks_lab = 0:0.1:1.0

# ============================================================
# ALPHAS
# ============================================================
alphas = [
    1e-12,
    1.6e-11,
    1e-9,
    1e-7
]

# ============================================================
# RNG THREADS
# ============================================================
rngs = [
    MersenneTwister(1234 + i)
    for i in 1:nthreads()
]

# ============================================================
# τ(x)
# ============================================================
function τ_exact(x, α)

    Nd = Ns * erfc(x / (2 * sqrt(D_Li * t_ann)))

    return pref_τ / (α * Nd)
end

function build_tau_grid(alpha)

    return [τ_exact(x, alpha) for x in xgrid]

end

@inline function τ_hole(x, τgrid)

    idx = Int(clamp(floor(x / LUT_dx) + 1, 1, length(τgrid)))

    return τgrid[idx]
end

# ============================================================
# trapping probability
# ============================================================
@inline function trapping_probability(x, τgrid)

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
    y = rand(rng) * Ly
    z = rand(rng) * Lz

    x <= dx && return (false, 0.0)

    @inbounds for i in 1:Nt

        x += step_sigma * randn(rng)
        y += step_sigma * randn(rng)
        z += step_sigma * randn(rng)

        # fast clamp
        if y < 0
            y = 0
        elseif y > Ly
            y = Ly
        end

        if z < 0
            z = 0
        elseif z > Lz
            z = Lz
        end

        # escape
        x <= 0.001 && return (false, 0.0)

        # trapping
        rand(rng) < trapping_probability(x, τgrid) &&
            return (false, 0.0)

        # collected
        x >= FCCD_cm &&
            return (true, i * Δt)
    end

    return (false, 0.0)
end

# ============================================================
# run depth
# ============================================================
function run_depth(x0, N, rng, τgrid)

    collected = 0

    times = Vector{Float64}(undef, N)
    n_times = 0

    @inbounds for i in 1:N

        ok, t = propagate_charge(x0, rng, τgrid)

        if ok

            collected += 1

            n_times += 1
            times[n_times] = t
        end
    end

    return collected, times[1:n_times]
end

# ============================================================
# MAIN PLOTS
# ============================================================
p1 = plot(
    xlabel="Depth (mm)",
    ylabel="CCE",
    frame=:box,
    xticks=xticks_lab,
    legend=:bottomright,
    lw=2
)

p2 = plot(
    xlabel="Depth (mm)",
    ylabel="t₉₀ - t₁₀ [μs]",
    frame=:box,
    xticks=xticks_lab,
    legend=:topleft,
    lw=2
)

# ============================================================
# GLOBAL PROGRESS
# ============================================================
total_steps = length(alphas) * Nx * N_repeat

pbar = Progress(total_steps)

# ============================================================
# LOOP OVER α
# ============================================================
for alpha_val in alphas

    println("\nRunning α = $alpha_val")

    τgrid = build_tau_grid(alpha_val)

    times_by_depth =
        Vector{Vector{Vector{Float64}}}(undef, Nx)

    N_matrix = zeros(Int, Nx, N_repeat)

    # ========================================================
    # PARALLEL SIMULATION
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
        fillalpha=0.5,
        label="α = $(alpha_val)",
    )

    # ========================================================
    # t90 - t10
    # ========================================================
    x_all = Float64[]
    dt_all = Float64[]

    for (i, x0) in enumerate(x_pos)

        reps = times_by_depth[i]

        for rep in reps

            if length(rep) > 5

                t10 = quantile(rep, 0.10)
                t90 = quantile(rep, 0.90)

                push!(x_all, x0 * 10)
                push!(dt_all, t90 - t10)
            end
        end
    end

    x_unique = unique(x_all)

    mean_dt = [
        mean(dt_all[x_all.==x] .* 0.001)
        for x in x_unique
    ]

    std_dt = [
        std(dt_all[x_all.==x] .* 0.001)
        for x in x_unique
    ]

    plot!(
        p2,
        x_unique,
        mean_dt,
        ribbon=std_dt,
        fillalpha=0.5,
        label="α = $(alpha_val)",
        legend=:topleft
    )
end

# ============================================================
# FINAL COMBINED FIGURE
# ============================================================
pfinal = plot(
    p1,
    p2,
    layout=(1, 2),
    size=(600, 450),
    dpi=300
)

# ============================================================
# SAVE
# ============================================================
mkpath("plot")

savefig(
    pfinal,
    "plot/CCE_and_t90t10_comparison.png"
)

println("\nALL DONE")
println("Saved:")
println("plot/CCE_and_t90t10_comparison.png")