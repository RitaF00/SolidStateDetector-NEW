using Random
using Plots
using SpecialFunctions
using ProgressMeter
using Statistics
using Base.Threads

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
#const α_global = 1.6e-11

const α_global = 0


const t_ann = 18 * 60
const T_ann = 623

const R = 1.98
const H = 11800
const D0 = 2.5e-3

const Ns = 10^(21.27 - 2610 / T_ann)
const D_Li = D0 * exp(-H / (R * T_ann))

const D = 28.9
const Δt = 1.0
const Nt = 40000  # 10 μs
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
# LOOKUP TABLE τ(x)
# ============================================================
const LUT_dx = step_sigma / 4
const xgrid = 0:LUT_dx:FCCD_cm

function τ_exact(x, α)

    Nd = Ns * erfc(x / (2 * sqrt(D_Li * t_ann)))

    return pref_τ / (α * Nd)

end

const τgrid = [τ_exact(x, α_global) for x in xgrid]

const step = 0.005

# ============================================================
# τ(x) via lookup
# ============================================================
@inline function τ_hole(x)

    idx = Int(clamp(floor(x / LUT_dx) + 1, 1, length(τgrid)))

    return τgrid[idx]
end

# ============================================================
# trapping probability
# ============================================================
@inline function trapping_probability(x)

    if α_global == 0.0
        return 0.0
    end

    if x >= x_RDR
        return 0.0
    end

    τ_ns = τ_hole(x) * 1e9

    if τ_ns <= 0
        return 1.0
    end

    return 1 - exp(-Δt / τ_ns)
end

# ============================================================
# propagazione singola carica
# ============================================================
@inline function propagate_charge(x0, rng)

    x = x0
    y = Ly / 2
    z = Lz / 2

    if x <= dx
        return (false, 0.0)
    end

    @inbounds for i in 1:Nt

        x += step_sigma * randn(rng)
        y += step_sigma * randn(rng)
        z += step_sigma * randn(rng)

        # clamp
        y = clamp(y, 0, Ly)
        z = clamp(z, 0, Lz)

        # exit p+
        if x <= step
            return (false, i * Δt)
        end

        # trapping
        if rand(rng) < trapping_probability(x)
            return (false, i * Δt)
        end

        # collection
        if x >= FCCD_cm
            return (true, i * Δt)
        end
    end

    return (false, Nt * Δt)
end
# ============================================================
# simulazione profondità
# ============================================================
function run_depth(x0, N, rng)

    collected = 0

    times = Vector{Float64}(undef, N)
    n_times = 0

    @inbounds for i in 1:N

        ok, t = propagate_charge(x0, rng)

        if ok
            collected += 1

            n_times += 1
            times[n_times] = t

        else
            n_times += 1
            times[n_times] = t

        end
    end

    return collected, times[1:n_times]
end

# ============================================================
# PARAMETRI SIMULAZIONE
# ============================================================

x_pos = step:step:FCCD_cm+step

N_charges = 100
N_repeat = 100

Nx = length(x_pos)

times_by_depth =
    Vector{Vector{Vector{Float64}}}(undef, Nx)

N_matrix = zeros(Int, Nx, N_repeat)

# ============================================================
# THREAD RNG
# ============================================================
rngs = [
    MersenneTwister(1234 + i)
    for i in 1:nthreads()
]

# ============================================================
# SIMULAZIONE PARALLELA
# ============================================================
pbar = Progress(Nx * N_repeat)

Threads.@threads for i in eachindex(x_pos)

    rng = rngs[threadid()]

    x0 = x_pos[i]

    rep_times =
        Vector{Vector{Float64}}(undef, N_repeat)

    for j in 1:N_repeat

        collected, times =
            run_depth(x0, N_charges, rng)

        N_matrix[i, j] = collected

        rep_times[j] = times

        next!(pbar)
    end

    times_by_depth[i] = rep_times
end

println("\nSIMULATION DONE")

# ============================================================
# SPREAD (t90 - t10)
# ============================================================
x_all = Float64[]
dt_all = Float64[]

t10_all = Float64[]
t90_all = Float64[]

for (i, x0) in enumerate(x_pos)

    reps = times_by_depth[i]

    for rep in reps

        if length(rep) > 5

            t10 = quantile(rep, 0.10)
            t90 = quantile(rep, 0.90)

            push!(x_all, x0 * 10)
            push!(dt_all, t90 - t10)
            push!(t10_all, t10)
            push!(t90_all, t90)
        end
    end
end

# ============================================================
# PLOT SPREAD
# ============================================================
blu_default = RGB(0.0, 0.6056, 0.9787)

xticks_lab = 0:0.1:1.4

p1 = scatter(
    x_all,
    dt_all .* 0.001,
    ms=3,
    alpha=0.6,
    color=:grey,
    xlabel="Depth (mm)",
    ylabel="t₉₀ - t₁₀ [μs]",
    label="",
    xticks=xticks_lab,
    frame=:box,
    size=(500, 400),
    dpi=300
)

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
    p1,
    x_unique,
    mean_dt,
    yerr=std_dt,
    lw=2,
    markercolor=blu_default,
    errcolor=blu_default,
    linecolor=blu_default,
    label="mean ± std"
)

# ============================================================
# CCE
# ============================================================
CCE = N_matrix ./ N_charges

CCE_mean = vec(mean(CCE, dims=2))
CCE_std = vec(std(CCE, dims=2))

p2 = plot(
    x_pos .* 10,
    CCE_mean,
    ribbon=CCE_std,
    fillalpha=0.2,
    lw=2,
    xlabel="Depth (mm)",
    ylabel="CCE",
    label="CCE",
    xticks=xticks_lab,
    frame=:box,
    size=(350, 500),
    dpi=300
)

# ============================================================
# DISTRIBUZIONI
# ============================================================
function adaptive_bins(data; min_bins=5, max_bins=200)
    n = length(data)
    n < 2 && return min_bins

    q25 = quantile(data, 0.25)
    q75 = quantile(data, 0.75)
    iqr = q75 - q25

    iqr == 0 && return min_bins

    h = 2 * iqr / n^(1 / 3)
    h <= 0 && return min_bins

    nbins = Int(clamp(ceil((maximum(data) - minimum(data)) / h), min_bins, max_bins))
    return nbins
end

n = length(x_pos)

ncols = 4
nrows = ceil(Int, n / ncols)

p3 = plot(
    layout=(nrows, ncols),
    size=(1200, 900),)

for i in 1:n
    all_data = Float64[]
    for v in times_by_depth[i]
        append!(all_data, v)
    end
    isempty(all_data) && continue

    data_μs = all_data ./ 1000

    nbins = adaptive_bins(data_μs)

    t10 = quantile(data_μs, 0.10)
    t90 = quantile(data_μs, 0.90)

    histogram!(
        p3,
        data_μs,
        bins=nbins,
        alpha=0.6,
        subplot=i,
        label="",
        frame=:box,
        xlim=(minimum(data_μs), t90 * 4),
        title="x = $(round(x_pos[i]*10, digits=3)) mm",
        xlabel="time (μs)"
    )

    vline!(p3, [t10], subplot=i, lc=:red, lw=2, ls=:dash, label="10%")
    vline!(p3, [t90], subplot=i, lc=:red, lw=2, ls=:dash, label="90%")
end


# ============================================================
# SAVE
# ============================================================
mkpath("plottini")

savefig(p1, "plottini/spread_scatter_precipitates_const.png")
savefig(p2, "plottini/CCE_precipitates_const.png")
savefig(p3, "plottini/tau_distributions_precipitates_const.png")




p4 = scatter(
    x_all,
    t10_all .* 0.001,
    ms=3,
    alpha=0.6,
    color=:grey,
    xlabel="Depth (mm)",
    ylabel=" t₁₀ [μs]",
    label="",
    xticks=xticks_lab,
    frame=:box,
    size=(500, 400),
    dpi=300
)

x_unique = unique(x_all)

mean_dt = [
    mean(t10_all[x_all.==x] .* 0.001)
    for x in x_unique
]

std_dt = [
    std(t10_all[x_all.==x] .* 0.001)
    for x in x_unique
]

plot!(
    p4,
    x_unique,
    mean_dt,
    yerr=std_dt,
    lw=2,
    markercolor=blu_default,
    errcolor=blu_default,
    linecolor=blu_default,
    label="mean ± std"
)


p5 = scatter(
    x_all,
    t90_all .* 0.001,
    ms=3,
    alpha=0.6,
    color=:grey,
    xlabel="Depth (mm)",
    ylabel=" t₉₀  [μs]",
    label="",
    xticks=xticks_lab,
    frame=:box,
    size=(500, 400),
    dpi=300
)

x_unique = unique(x_all)

mean_dt = [
    mean(t90_all[x_all.==x] .* 0.001)
    for x in x_unique
]

std_dt = [
    std(t90_all[x_all.==x] .* 0.001)
    for x in x_unique
]

plot!(
    p5,
    x_unique,
    mean_dt,
    yerr=std_dt,
    lw=2,
    markercolor=blu_default,
    errcolor=blu_default,
    linecolor=blu_default,
    label="mean ± std"
)

p = plot(p4, p5, layout=(1, 2), size=(900, 400), legend=false)
savefig(p, "t10+t90.png")

println("\nALL DONE")
println("Plots saved in ./plottini/")