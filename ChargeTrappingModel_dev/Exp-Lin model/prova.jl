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
const α_global = 8e-9

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

const step_sigma = sqrt(2 * D * Δt) * 1e-4

# ============================================================
# LOOKUP τ(x)
# ============================================================
const LUT_dx = step_sigma / 4
const xgrid = 0:LUT_dx:FCCD_cm

function τ_exact(x, α)
    Nd = Ns * erfc(x / (2 * sqrt(D_Li * t_ann)))
    return pref_τ / (α * Nd)
end

const τgrid = [τ_exact(x, α_global) for x in xgrid]

@inline function τ_hole(x)
    idx = Int(clamp(floor(x / LUT_dx) + 1, 1, length(τgrid)))
    return τgrid[idx]
end

# ============================================================
# TRAPPING
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
# PROPAGAZIONE
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

        y = clamp(y, 0, Ly)
        z = clamp(z, 0, Lz)

        if x <= 0
            return (false, i * Δt)
        end

        if rand(rng) < trapping_probability(x)
            return (false, i * Δt)
        end

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

    for i in 1:N
        ok, t = propagate_charge(x0, rng)
        times[i] = t
        ok && (collected += 1)
    end

    return collected, times
end

# ============================================================
# SIMULAZIONE
# ============================================================
x_pos = 0:0.0025:(FCCD_cm+5*0.0025)

N_charges = 100
N_repeat = 100
Nx = length(x_pos)

N_matrix = zeros(Int, Nx, N_repeat)

rngs = [MersenneTwister(1234 + i) for i in 1:nthreads()]

cache_file = "CCE_cache.jld2"
force_recompute = false

if isfile(cache_file) && !force_recompute
    @load cache_file N_matrix
else
    pbar = Progress(Nx * N_repeat)

    Threads.@threads for i in eachindex(x_pos)
        rng = rngs[threadid()]
        x0 = x_pos[i]

        for j in 1:N_repeat
            collected, _ = run_depth(x0, N_charges, rng)
            N_matrix[i, j] = collected
            next!(pbar)
        end
    end

    @save cache_file N_matrix
end

# ============================================================
# CCE
# ============================================================
CCE = N_matrix ./ N_charges
CCE_mean = vec(mean(CCE, dims=2))
CCE_std = vec(std(CCE, dims=2))
CCE_std[CCE_std.==0.0] .= 1e-6

x_mm = x_pos .* 1

# ============================================================
# FIT REGION (robusto)
# ============================================================
eps_CCE = 1e-4
idx0 = findfirst(x -> x > eps_CCE, CCE_mean)
x_start_exp = isnothing(idx0) ? 1 : idx0

mask_exp = (1:length(x_mm) .>= x_start_exp) .& (x_mm .<= x_RDR)
mask_lin = x_mm .>= x_RDR

x_exp = x_mm[mask_exp]
y_exp = CCE_mean[mask_exp]

x_lin = x_mm[mask_lin]
y_lin = CCE_mean[mask_lin]

x_DL = x_mm[x_start_exp]

# ============================================================
# MODELLO ESPONENZIALE CON x0
# ============================================================
CCE_exp(x, p) = p[1] .* (exp.(p[2] .* (x .- p[5]) .^ 2 .+
                              p[3] .* (x .- p[5]) .+
                              p[4]) .- 1)

p0exp = [1e-6, 10.0, -10.0, 0.0, x_DL]
fit_exp = curve_fit(CCE_exp, x_exp, y_exp, p0exp)
param_exp = fit_exp.param

# ============================================================
# DERIVATA (C¹)
# ============================================================
function dCCE_exp(x, p)
    A, B, C, D, x0 = p

    u = x - x0
    g = B * u^2 + C * u + D
    dg = 2B * u + C

    return A * exp(g) * dg
end

# ============================================================
# MATCH C¹
# ============================================================
y_match = CCE_exp([x_RDR], param_exp)[1]
m_match = dCCE_exp(x_RDR, param_exp)

# ============================================================
# RETTA C¹
# ============================================================
CCE_lin(x) = m_match .* (x .- x_RDR) .+ y_match

# ============================================================
# OUTPUT
# ============================================================
println("===================================")
println("FIT ESPONENZIALE CON x0")
println("A,B,C,D,x0 = ", param_exp)
println("===================================")
println("C¹ MATCH")
println("y_match = ", y_match)
println("slope = ", m_match)

# ============================================================
# PLOT
# ============================================================
x_smooth_exp = range(minimum(x_exp), maximum(x_exp), length=300)
x_smooth_lin = range(x_RDR, maximum(x_lin), length=300)

p_fit = plot(x_mm .* 10, CCE_mean, ribbon=CCE_std, label="CCE")
plot!(p_fit, x_smooth_exp .* 10, CCE_exp(x_smooth_exp, param_exp),
    ls=:dash, label="exp fit")
plot!(p_fit, x_smooth_lin .* 10, CCE_lin(x_smooth_lin),
    ls=:dashdot, label="C¹ continuation")

scatter!(p_fit, [x_RDR * 10], [y_match], label="C¹ match")

display(p_fit)