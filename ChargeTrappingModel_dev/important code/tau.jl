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
const α_global = 1.6e-11

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
# TRAPPING PROBABILITY
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
# PROPAGAZIONE CARICA
# ============================================================
@inline function propagate_charge(x0, rng)

    x = x0
    y = Ly / 2
    z = Lz / 2

    if x <= dx
        return (false, 0.0)
    end

    for i in 1:Nt

        x += step_sigma * randn(rng)
        y += step_sigma * randn(rng)
        z += step_sigma * randn(rng)

        y = clamp(y, 0, Ly)
        z = clamp(z, 0, Lz)

        if x <= step_sigma
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
# SIMULAZIONE PER PROFONDITÀ
# ============================================================
function run_depth(x0, N, rng)

    collected = 0
    times = Float64[]

    for i in 1:N
        ok, t = propagate_charge(x0, rng)
        collected += ok ? 1 : 0
        push!(times, t)
    end

    return collected, times
end

# ============================================================
# PARAMETRI SIMULAZIONE CCE
# ============================================================
x_pos = step_sigma:0.005:FCCD_cm
N_charges = 100
N_repeat = 100

Nx = length(x_pos)

N_matrix = zeros(Int, Nx, N_repeat)

rng = MersenneTwister(1234)

# ============================================================
# SIMULAZIONE CCE
# ============================================================
pbar = Progress(Nx * N_repeat)

for (i, x0) in enumerate(x_pos)
    for j in 1:N_repeat

        collected, _ = run_depth(x0, N_charges, rng)
        N_matrix[i, j] = collected

        next!(pbar)
    end
end

println("\nSIMULATION DONE")

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
    frame=:box,
    size=(400, 500),
    dpi=300
)

display(p2)