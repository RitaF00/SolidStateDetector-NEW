using Unitful
using Plots
using SpecialFunctions

"""
In questo codice vedo come si comporta τ_hole(x) al variare di α, e come questo influenza la probabilità di trapping.
dipende solo da α, Nd(x) e v_th, che a sua volta dipende da T_diff. Quindi vedo come varia τ_hole(x) al variare di α, e come questo influenza la probabilità di trapping.
NON dipende dalla diffusione.
Valori corretti da mettere nel codice con CCE
"""
# =========================
# ⚙️ PARAMETRI
# =========================
T = Float64

det_r = 10u"mm"
pn_r = 8.957282u"mm"


# =========================
# ⚙️ PARAMETRI
# =========================
T_diff = 90u"K"
Δt = 1u"ns"

t_ann = 18 * 60u"s"
T_ann = 623u"K"

# =========================
# 📌 COSTANTI
# =========================
m0 = 9.11e-31u"kg"
kB = 1.380649e-23u"J/K"

# Litio
r_Li = 0.002u"cm"

R = 1.98u"cal/(mol*K)"
H = 11800u"cal/mol"
D0 = 2.5e-3u"cm^2/s"

# Concentrazione e diffusione
Ns = 10^(21.27 - 2610 / ustrip(T_ann)) * u"cm^-3"
D_Li = D0 * exp(-H / (R * T_ann))

# =========================
# ⏳ τ_hole(x)
# =========================
function τ_hole(x, α)
    m_eff = 0.21 * m0

    # velocità termica (m/s → convertita in cm/s)
    v_th = sqrt(3 * kB * T_diff / m_eff)
    v_th = uconvert(u"cm/s", v_th)

    σ = π * r_Li^2

    Nd = Ns .* erfc.(x ./ (2 * sqrt(D_Li * t_ann)))
    Nd = uconvert.(u"cm^-3", Nd)

    τ = 1.0 ./ (α .* Nd .* v_th .* σ)

    return uconvert.(u"s", τ)
end

# =========================
# 🎯 Probabilità di trapping
# =========================
function P_trapping(Δt, x, α)
    τ = τ_hole(x, α)
    return 1 .- exp.(-Δt ./ τ)
end

# =========================
# 📈 RANGE
# =========================
x = (0:0.001:0.11) .* u"cm"
alphas = [1e-8, 1e-9, 1e-10, 1e-11, 1e-12]

colors = [:deepskyblue, :royalblue1, :mediumpurple2, :deeppink, :orange]

# =========================
# 📊 PLOT τ(x)
# =========================
p = plot(
    ustrip.(u"mm", x),
    ustrip.(u"ns", τ_hole(x, alphas[1])),
    lw=3,
    color=colors[1],
    label="α = $(alphas[1])",
    yscale=:log10,
    frame=:box,
    xlabel="depth [mm]",
    ylabel="τ [ns]",
    legend=:topleft,
    size=(600, 400)
)

for i in 2:length(alphas)
    plot!(
        p,
        ustrip.(u"mm", x),
        ustrip.(u"ns", τ_hole(x, alphas[i])),
        lw=3,
        color=colors[i],
        label="α = $(alphas[i])"
    )
end

# linea τ costante = 800 ns
plot!(
    p,
    ustrip.(u"mm", x),
    fill(800, length(x)),
    linestyle=:dash,
    color=:red,
    lw=2,
    label="τ = 800 ns"
)

display(p)

p2 = plot(
    ustrip.(u"mm", x),
    P_trapping(Δt, x, alphas[1]),
    lw=3,
    color=colors[1],
    label="α = $(alphas[1])",
    frame=:box,
    xlabel="depth [mm]",
    ylabel="P trapping",
    legend=:bottomright,
    size=(600, 400)
)


for i in 2:length(alphas)
    plot!(
        p2,
        ustrip.(u"mm", x),
        P_trapping(Δt, x, alphas[i]),
        lw=3,
        color=colors[i],
        label="α = $(alphas[i])",
        yscale=:log10
    )
end

display(p2)


function Nd_profile(x)
    arg = ustrip.(x ./ (2 * sqrt(D_Li * t_ann)))
    Nd = Ns .* erfc.(arg)
    return uconvert.(u"cm^-3", Nd)
end

Nd = Nd_profile(x)

p3 = plot(
    ustrip.(u"mm", x),
    ustrip.(Nd .* alphas[1]),
    lw=3,
    color=colors[1],
    label="α = $(alphas[1])",
    frame=:box,
    xlabel="depth [mm]",
    ylabel="α · Nd [cm⁻³]",
    yscale=:log10,
    legend=:topright,
    size=(600, 400)
)

for i in 2:length(alphas)
    plot!(
        p3,
        ustrip.(u"mm", x),
        ustrip.(Nd .* alphas[i]),
        lw=3,
        color=colors[i],
        label="α = $(alphas[i])"
    )
end

display(p3)