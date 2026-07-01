using Random
using Plots
using SpecialFunctions
using LinearAlgebra
using ProgressMeter
using Statistics
using Base.Threads
using JLD2
using LsqFit

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
const α_global = 8e-13

const t_ann = 18 * 60
const T_ann = 623

const R = 1.98
const H = 11800
const D0 = 2.5e-3

const Ns = 10^(21.27 - 2610 / T_ann)
const D_Li = D0 * exp(-H / (R * T_ann))

const D = 28.9
const Δt = 1.0
const Nt = 40000
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

# step spaziale random walk
const step_sigma = sqrt(2 * D * Δt) * 1e-4

# ============================================================
# LOOKUP TABLE τ(x)
# ============================================================
const LUT_dx = step_sigma / 4
const xgrid = 0:LUT_dx:FCCD_cm

function τ_exact(x, α)

    Nd = Ns * erfc(
        x / (2 * sqrt(D_Li * t_ann))
    )

    return pref_τ / (α * Nd)
end

const τgrid = [
    τ_exact(x, α_global)
    for x in xgrid
]

# soglia raccolta
const step = 0.0025

@inline function τ_hole(x)

    idx = Int(
        clamp(
            floor(x / LUT_dx) + 1,
            1,
            length(τgrid)
        )
    )

    return τgrid[idx]
end

# ============================================================
# TRAPPING PROBABILITY
# ============================================================
@inline function trapping_probability(x)

    if α_global == 0.0 || x >= x_RDR
        return 0.0
    end

    τ_ns = τ_hole(x) * 1e9

    τ_ns <= 0 && return 1.0

    return 1 - exp(-Δt / τ_ns)
end

# ============================================================
# PROPAGAZIONE SINGOLA CARICA
# ============================================================
@inline function propagate_charge(x0, rng)

    x = x0
    y = Ly / 2
    z = Lz / 2

    # morto immediatamente
    if x <= dx
        return (false, 0.0)
    end

    @inbounds for i in 1:Nt

        # random walk 3D
        x += step_sigma * randn(rng)
        y += step_sigma * randn(rng)
        z += step_sigma * randn(rng)

        # riflessione laterale
        y = clamp(y, 0, Ly)
        z = clamp(z, 0, Lz)

        # perso
        if x <= step
            return (false, i * Δt)
        end

        # trapping
        if rand(rng) < trapping_probability(x)
            return (false, i * Δt)
        end

        # raccolto
        if x >= FCCD_cm
            return (true, i * Δt)
        end
    end

    return (false, Nt * Δt)
end

# ============================================================
# SIMULAZIONE PROFONDITÀ
# ============================================================
function run_depth(x0, N, rng)

    collected = 0

    times = Vector{Float64}(undef, N)
    n_times = 0

    @inbounds for i in 1:N

        ok, t = propagate_charge(x0, rng)

        n_times += 1
        times[n_times] = t

        ok && (collected += 1)
    end

    return collected, times[1:n_times]
end

# ============================================================
# PARAMETRI SIMULAZIONE
# ============================================================
x_pos = 0:step:(FCCD_cm+5*step)

N_charges = 100
N_repeat = 100

Nx = length(x_pos)

times_by_depth =
    Vector{Vector{Vector{Float64}}}(undef, Nx)

N_matrix = zeros(Int, Nx, N_repeat)

rngs = [
    MersenneTwister(1234 + i)
    for i in 1:nthreads()
]

# ============================================================
# CACHE CCE
# ============================================================
cache_file = "CCE_cache.jld2"

force_recompute_CCE = false

if isfile(cache_file) && !force_recompute_CCE

    @load cache_file N_matrix times_by_depth

else

    pbar = Progress(Nx * N_repeat)

    Threads.@threads for i in eachindex(x_pos)

        rng = rngs[threadid()]

        x0 = x_pos[i]

        rep_times =
            Vector{Vector{Float64}}(
                undef,
                N_repeat
            )

        for j in 1:N_repeat

            collected, times =
                run_depth(
                    x0,
                    N_charges,
                    rng
                )

            N_matrix[i, j] = collected

            rep_times[j] = times

            next!(pbar)
        end

        times_by_depth[i] = rep_times
    end

    @save cache_file N_matrix times_by_depth
end

# ============================================================
# CCE
# ============================================================
CCE = N_matrix ./ N_charges

CCE_mean = vec(mean(CCE, dims=2))
CCE_std = vec(std(CCE, dims=2))

# evita divisioni per zero
CCE_std[CCE_std.==0.0] .= 1e-6

# ============================================================
# REGIONI FIT
# ============================================================
x_mm = x_pos .* 1

mask_exp =
    (x_mm .>= 0) .&
    (x_mm .<= x_RDR)

mask_lin =
    (x_mm .>= x_RDR) .&
    (x_mm .<= FCCD_cm)

x_exp = x_mm[mask_exp]
y_exp = CCE_mean[mask_exp]

x_lin = x_mm[mask_lin]
y_lin = CCE_mean[mask_lin]

# ============================================================
# MODELLI
# ============================================================

# passa per (0,0)
CCE_exp(x, p) =
    p[1] .* (
        exp.(p[2] .* x .* x .+ p[3] .* x) .- 1
    )

# ============================================================
# FIT ESPONENZIALE
# ============================================================
p0exp = [1e-6, 100.0, 100.0]

fit_exp =
    curve_fit(
        CCE_exp,
        x_exp,
        y_exp,
        p0exp
    )

param_exp = fit_exp.param
se_exp = stderror(fit_exp)

# ============================================================
# MATCH CONTINUO NEL PUNTO x_RDR
# ============================================================
y_match =
    CCE_exp(
        [x_RDR],
        param_exp
    )[1]

# ============================================================
# MODELLO LINEARE VINCOLATO
# impone continuità:
#
# CCE_lin(x_RDR) = CCE_exp(x_RDR)
# ============================================================
CCE_lin(x, p) =
    p[1] .* (x .- x_RDR) .+ y_match

# unico parametro libero = slope
p0lin = [1.0]

# ============================================================
# FIT LINEARE
# ============================================================
fit_lin =
    curve_fit(
        CCE_lin,
        x_lin,
        y_lin,
        p0lin
    )

param_lin = fit_lin.param
se_lin = stderror(fit_lin)

# ============================================================
# STAMPA RISULTATI
# ============================================================
println()
println("================================================")
println("EXP FIT")
println("================================================")
println("params = ", param_exp)
println("errors = ", se_exp)

println()
println("================================================")
println("LIN FIT (continuo in x_RDR)")
println("================================================")
println("params = ", param_lin)
println("errors = ", se_lin)

println()
println("y_match = ", y_match)

# ============================================================
# CURVE LISCE
# ============================================================
x_smooth_exp =
    range(
        minimum(x_exp),
        maximum(x_exp),
        length=300
    )

x_smooth_lin =
    range(
        minimum(x_lin),
        maximum(x_lin),
        length=300
    )

y_exp_smooth =
    CCE_exp(
        x_smooth_exp,
        param_exp
    )

y_lin_smooth =
    CCE_lin(
        x_smooth_lin,
        param_lin
    )

# ============================================================
# RESIDUI
# ============================================================
y_exp_fit =
    CCE_exp(
        x_exp,
        param_exp
    )

y_lin_fit =
    CCE_lin(
        x_lin,
        param_lin
    )

σ_exp = CCE_std[mask_exp]
σ_lin = CCE_std[mask_lin]

res_exp =
    (y_exp .- y_exp_fit) ./ σ_exp

res_lin =
    (y_lin .- y_lin_fit) ./ σ_lin

residuals =
    vcat(
        res_exp,
        res_lin
    )

x_res =
    vcat(
        x_exp,
        x_lin
    )

# ============================================================
# PLOT CCE
# ============================================================
p_fit = plot(
    x_mm .* 10,
    CCE_mean,
    ribbon=CCE_std,
    fillalpha=0.2,
    lw=2,
    xlabel="Depth (mm)",
    ylabel="CCE",
    label="CCE",
    frame=:box,
    dpi=300,
    size=(450, 400)
)

# exp fit
plot!(
    p_fit,
    x_smooth_exp .* 10,
    y_exp_smooth,
    lw=3,
    ls=:dash,
    color=:orange,
    label="Exp fit"
)

# linear fit
plot!(
    p_fit,
    x_smooth_lin .* 10,
    y_lin_smooth,
    lw=3,
    ls=:dashdot,
    color=:red,
    label="Linear fit"
)

# punto di giunzione
scatter!(
    p_fit,
    [x_RDR * 10],
    [y_match],
    ms=5,
    color=:black,
    label="Match point"
)

# ============================================================
# PLOT RESIDUI
# ============================================================
p_res = scatter(
    x_res .* 10,
    residuals,
    ms=2,
    alpha=0.8,
    color=:black,
    xlabel="Depth (mm)",
    ylabel="Normalized residuals (σ)",
    frame=:box,
    label=""
)

hline!(
    p_res,
    [0.0],
    lc=:black,
    lw=1,
    ls=:dash,
    label=""
)

plot!(
    p_res,
    x_res .* 10,
    fill(1, length(x_res)),
    fillrange=fill(-1, length(x_res)),
    fillalpha=0.25,
    linecolor=:transparent,
    fillcolor=:green,
    label=""
)

plot!(
    p_res,
    x_res .* 10,
    fill(2, length(x_res)),
    fillrange=fill(-2, length(x_res)),
    fillalpha=0.18,
    linecolor=:transparent,
    fillcolor=:orange,
    label=""
)

plot!(
    p_res,
    x_res .* 10,
    fill(3, length(x_res)),
    fillrange=fill(-3, length(x_res)),
    fillalpha=0.12,
    linecolor=:transparent,
    fillcolor=:red,
    label=""
)

scatter!(
    p_res,
    x_res .* 10,
    residuals,
    ms=2,
    alpha=0.8,
    color=:black,
    label=""
)

# ============================================================
# PLOT FINALE
# ============================================================
p_final = plot(
    p_fit,
    p_res,
    layout=@layout([a{0.7h}; b{0.3h}]),
    size=(400, 900)
)

#display(p_final)

# ============================================================
# SAVE
# ============================================================
mkpath("plottini")

savefig(
    p_final,
    "plottini/CCE_precipitates_const.png"
)



#=


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
# SAVE
# ============================================================
mkpath("plottini")
savefig(p2, "plottini/CCE_precipitates_const.png")
=#