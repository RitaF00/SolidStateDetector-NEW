using Random
using Plots
using SpecialFunctions
using ProgressMeter
using Distributions
using JSON
using Statistics
using LaTeXStrings



T_diff = 90 #K
r_Li = 0.002  # cm, 20 µm
Δt = 1 # ns

t_ann = 18 * 60  # s
T_ann = 623  # K

R = 1.98
H = 11800
D0 = 2.5e-3

Ns = 10^(21.27 - 2610 / T_ann)   #cm^-3
D_Li = D0 * exp(-H / (R * T_ann))


function τ_hole(x, α)
    m0 = 9.11e-31 #kg
    m_eff_hole = 0.21 * m0
    kB = 1.38e-23 #J/K
    v_th = [(3 * kB * T_diff / m_eff_hole)^0.5] .* 1e-2  # cm/s
    σ_trap = π * (r_Li)^2
    Nd_val = Ns * erfc.(x ./ (2 * sqrt(D_Li * t_ann)))   #cm^.3

    τ_h = 1 ./ (α .* Nd_val .* v_th .* σ_trap)
    return τ_h
end


function P_trapping(Δt, x, α)
    # Δt --> ns
    τ_hole_ns = τ_hole(x, α) .* 1e9
    return 1 .- exp.(-Δt ./ τ_hole_ns)
end


x = 0:0.001:0.11 # cm
alphas = [1e-8, 1e-9, 1e-10, 1e-11, 1e-12]
alphas = [1.6e-11]
colors = [:deepskyblue, :royalblue1, :mediumpurple2, :deeppink, :orange]
exponents = -8:1
yticks = 10.0 .^ exponents
ylabels = [L"10^{%$e}" for e in exponents]
xticks = 0:0.1:1.1   # ogni 0.1 mm

#--- τ
exponents = 2:11
yticks = 10.0 .^ exponents
ylabels = [L"10^{%$e}" for e in exponents]
xticks = 0:0.1:1.1   # ogni 0.1 mm

p = plot(
    x .* 10,
    τ_hole(x, alphas[1]) .* 1e9,
    color=colors[1],
    lw=3,
    label="α = $(alphas[1])",
    yticks=(yticks, ylabels),
    xticks=xticks,
    yscale=:log10,
    frame=:box,
    xlabel="depth [mm]",
    ylabel="τ [ns]",
    legend=:topleft,
    legendfontsize=10
)

for i in 2:length(alphas)
    plot!(
        p,
        x .* 10,
        lw=3,
        τ_hole(x, alphas[i]) .* 1e9,
        color=colors[i],
        label="α = $(alphas[i])"
    )
end

plot!(x .* 10, fill((800), length(x)),
    linestyle=:dash,
    color=:red,
    lw=2,
    label="τ = 800 ns")

display(p)
# --- probabilità
#=
p = plot(
    x .* 10,
    P_trapping(Δt, x, alphas[1]) .* 100,
    color=colors[1],
    lw=3,
    label="α = $(alphas[1])",
    yticks=(yticks, ylabels),
    xticks=xticks,
    yscale=:log10,
    frame=:box,
    xlabel="depth [mm]",
    ylabel="Probability of trapping [%]",
    legend=:bottomleft,
    legendfontsize=10
)

for i in 2:length(alphas)
    plot!(
        p,
        x .* 10,
        lw=3,
        P_trapping(Δt, x, alphas[i]) .* 1e2,
        color=colors[i],
        label="α = $(alphas[i])"
    )
end

plot!(x .* 10, fill(1 - exp(-Δt / 800), length(x)),
    linestyle=:dash,
    color=:red,
    lw=2,
    label="τ = 800 ns")


display(p)
=#