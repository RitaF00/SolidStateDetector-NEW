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
alpha_list = [0, 1e-13, 1e-12, 1.6e-11, 5e-10, 1e-8,]

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

const LUT_dx = 1e-4
const xgrid = 0:LUT_dx:FCCD_cm

# ============================================================
# STORAGE RISULTATI GLOBALI
# ============================================================
CCE_store = Dict{Float64,Any}()
spread_store = Dict{Float64,Any}()

# ============================================================
# SIMULAZIONE PER OGNI α
# ============================================================

times_store = Dict{Float64,Any}()
t10_store = Dict{Float64,Any}()
t90_store = Dict{Float64,Any}()


for α_global in alpha_list

    println("\n==============================")
    println("Running α = $α_global")
    println("==============================\n")

    # ========================================================
    # LOOKUP τ(x)
    # ========================================================


    τ_exact(x, α) = pref_τ / (α * (Ns * erfc(x / (2 * sqrt(D_Li * t_ann)))))
    #τ_exact = 1e-6
    τgrid = [τ_exact(x, α_global) for x in xgrid]
    #τgrid = [τ_exact for x in xgrid]

    @inline τ_hole(x) = τgrid[Int(clamp(floor(x / LUT_dx) + 1, 1, length(τgrid)))]


    @inline function trapping_probability(x)
        if x >= x_RDR
            return 0.0
        end

        τ_ns = τ_hole(x) * 1e9
        if τ_ns <= 0
            return 1.0
        end

        return 1 - exp(-Δt / τ_ns)
    end
    """

    @inline function trapping_probability(x)
        τ_ns = τ_hole(x) * 1e9
        if τ_ns <= 0
            return 1.0
        end

        return 1 - exp(-Δt / τ_ns)
    end
       """


    # ========================================================
    # TRASPORTO
    # ========================================================
    @inline function propagate_charge(x0, rng)

        x, y, z = x0, Ly / 2, Lz / 2

        if x <= dx
            return (false, 0.0)
        end

        @inbounds for i in 1:Nt

            x += step_sigma * randn(rng)
            y += step_sigma * randn(rng)
            z += step_sigma * randn(rng)

            y = clamp(y, 0, Ly)
            z = clamp(z, 0, Lz)

            if x <= 0.001
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

    function run_depth(x0, N, rng)

        collected = 0
        times = Float64[]

        for i in 1:N
            ok, t = propagate_charge(x0, rng)
            if ok
                collected += 1
            end
            push!(times, t)
        end

        return collected, times
    end

    # ========================================================
    # SIMULATION GRID
    # ========================================================
    step = 0.005
    x_pos = step:step:FCCD_cm-step

    N_charges = 100
    N_repeat = 500

    Nx = length(x_pos)

    times_by_depth = Vector{Vector{Vector{Float64}}}(undef, Nx)
    times_store[α_global] = times_by_depth
    N_matrix = zeros(Int, Nx, N_repeat)

    rngs = [MersenneTwister(1234 + i) for i in 1:nthreads()]

    pbar = Progress(Nx * N_repeat)

    Threads.@threads for i in eachindex(x_pos)

        rng = rngs[threadid()]
        x0 = x_pos[i]

        reps = Vector{Vector{Float64}}(undef, N_repeat)

        for j in 1:N_repeat

            collected, times = run_depth(x0, N_charges, rng)

            N_matrix[i, j] = collected
            reps[j] = times

            next!(pbar)
        end

        times_by_depth[i] = reps
    end

    println("\nDONE α = $α_global")

    # ========================================================
    # SPREAD
    # ========================================================
    x_all = Float64[]
    dt_all = Float64[]
    t10_map = Float64[]
    t90_map = Float64[]

    for (i, x0) in enumerate(x_pos)
        for rep in times_by_depth[i]
            if length(rep) > 5
                t10 = quantile(rep, 0.10)
                t90 = quantile(rep, 0.90)

                push!(x_all, x0 * 10)
                push!(dt_all, t90 - t10)
                push!(t10_map, t10)
                push!(t90_map, t90)
            end
        end
    end

    x_unique = unique(x_all)

    mean_dt = [mean(dt_all[x_all.==x] .* 0.001) for x in x_unique]
    std_dt = [std(dt_all[x_all.==x] .* 0.001) for x in x_unique]

    # ========================================================
    # CCE
    # ========================================================
    CCE = N_matrix ./ N_charges
    CCE_mean = vec(mean(CCE, dims=2))
    CCE_std = vec(std(CCE, dims=2))

    # ========================================================
    # STORE
    # ========================================================
    CCE_store[α_global] = (x_pos .* 10, CCE_mean, CCE_std)
    spread_store[α_global] = (x_unique, mean_dt, std_dt)
end

println("\nSIMULATIONS COMPLETE — plotting...")

# ============================================================
# PLOT CCE (OVERLAY)
# ============================================================
xticks_label = 0:0.1:1.1
p_CCE = plot(
    xlabel="Depth (mm)",
    ylabel="CCE",
    xticks=xticks_label,
    frame=:box,
    size=(650, 450),
    legend=:topright,
    dpi=300
)

for α in alpha_list
    x, y, σ = CCE_store[α]

    plot!(p_CCE, x, y, ribbon=σ, fillalpha=0.15, lw=2, label="α = $α")
end

# ============================================================
# PLOT SPREAD (OVERLAY)
# ============================================================
p_spread = plot(
    xlabel="Depth (mm)",
    ylabel="t₉₀ - t₁₀ [μs]",
    xticks=xticks_label,
    frame=:box,
    size=(650, 450),
    legend=:topright,
    dpi=300
)

for α in alpha_list
    x, m, s = spread_store[α]

    plot!(p_spread, x, m, ribbon=s, fillalpha=0.15, lw=2, label="α = $α")
end
x_pos = 0.005:0.005:FCCD_cm-0.005
p_t10 = plot(
    xlabel="Depth (mm)",
    ylabel="t₁₀ [μs]",
    xticks=xticks_label,
    frame=:box,
    size=(650, 450),
    legend=:topright,
    dpi=300
)

p_t90 = plot(
    xlabel="Depth (mm)",
    ylabel="t₉₀ [μs]",
    xticks=xticks_label,
    frame=:box,
    size=(650, 450),
    legend=:topright,
    dpi=300
)

for α in alpha_list

    times_by_depth = times_store[α]

    t10_vals = Float64[]
    t90_vals = Float64[]

    for i in 1:length(x_pos)

        reps = times_by_depth[i]

        vals10 = [quantile(rep, 0.10) for rep in reps if length(rep) > 5]
        vals90 = [quantile(rep, 0.90) for rep in reps if length(rep) > 5]

        push!(t10_vals, mean(vals10) * 0.001)
        push!(t90_vals, mean(vals90) * 0.001)
    end

    plot!(p_t10, x_pos .* 10, t10_vals, lw=2, label="α = $α")
    plot!(p_t90, x_pos .* 10, t90_vals, lw=2, label="α = $α")
end

# ============================================================
# SAVE
# ============================================================
mkpath("plottini")

savefig(p_CCE, "plottini/CCE_all_alpha.png")
savefig(p_spread, "plottini/spread_all_alpha.png")
savefig(p_t10, "plottini/t10_all_alpha.png")
savefig(p_t90, "plottini/t90_all_alpha.png")


mkpath("plottini/time_distributions")

for α in alpha_list

    println("Plotting distributions for α = $α")

    times_by_depth = times_store[α]

    n = length(x_pos)

    ncols = 4
    nrows = ceil(Int, n / ncols)

    p_dist = plot(
        layout=(nrows, ncols),
        size=(1200, 900),
        dpi=300
    )

    for i in 1:n

        all_data = Float64[]

        for v in times_by_depth[i]
            append!(all_data, v)
        end

        if isempty(all_data)
            continue
        end

        # ns → μs
        data_us = all_data ./ 1000

        t10 = quantile(data_us, 0.10)
        t90 = quantile(data_us, 0.90)

        histogram!(
            p_dist,
            data_us,
            bins=100,
            alpha=0.6,
            subplot=i,
            label="",
            frame=:box,
            xlim=(minimum(data_us), maximum(data_us)),
            title="x = $(round(x_pos[i] * 10, digits=3)) mm",
            xlabel="time (μs)"
        )

        vline!(
            p_dist,
            [t10],
            subplot=i,
            lc=:red,
            lw=2,
            ls=:dash,
            label="10%"
        )

        vline!(
            p_dist,
            [t90],
            subplot=i,
            lc=:blue,
            lw=2,
            ls=:dash,
            label="90%"
        )
    end

    # nome file pulito
    alpha_str = replace(string(α), "." => "p")

    savefig(
        p_dist,
        "plottini/time_distributions/time_distrib_alpha_$(alpha_str).png"
    )
end


println("DONE — plots saved in plottini/")