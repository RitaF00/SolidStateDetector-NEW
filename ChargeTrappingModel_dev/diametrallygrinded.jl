using CSV
using DataFrames
using Plots
using Measures

# =======================
# LOAD DATA
# =======================
function load_data(file)

    df = CSV.read(file, DataFrame; decimal=',')

    # colonne
    diametro = df[!, "diameter (cm)"]
    res_d = df[!, "resistance d (omh.cm)"]
    lunghezza = df[!, "length (cm)"]
    res_l = df[!, "resistance l (omh.cm)"]

    # conducibilità
    σ_l = skipmissing(1 ./ res_l) |> collect
    σ_d = skipmissing(1 ./ res_d) |> collect

    return lunghezza, diametro, σ_l, σ_d
end

# =======================
# FILE (3 DETECTOR V)
# =======================
files = ["V05612B.csv", "V01406A.csv", "P00748B.csv"]

data = load_data.(files)

L = [d[1] for d in data]
R = [d[2] for d in data]
σ = [d[3] for d in data]
σd = [d[4] for d in data]

# =======================
# PROFONDITÀ
# =======================
depth = [
    85 .- L[1],
    60.6 .- L[2],
    52.3 .- L[3]
]

# =======================
# COSTANTI FISICHE
# =======================
e = 1.602e-19
mu_e = 3.9e3
mu_h = 1.9e3
k = 8.617e-5
T = 300
Eg = 0.66

Nc = 1.04e19
Nv = 6.0e18

ni = sqrt(Nc * Nv) * exp(-Eg / (2 * k * T))

println(ni)
# =======================
# Nd FROM σ
# =======================
function Nd_from_sigma(σ)
    A = 4 * mu_e * mu_h * ni^2 * e^2

    rad = σ .^ 2 .- A
    rad = max.(rad, 0.0)

    term = sqrt.(rad)

    return (-σ .* (mu_e - mu_h) .+ (mu_e + mu_h) .* term) ./
           (2 * mu_e * mu_h * e)
end

Nd = Nd_from_sigma.(σ)


# =======================
# OUTPUT
# =======================
using Printf

for i in 1:3
    val = maximum(Nd[i])
    @printf "Saturation Nd for %s : %.2e cm^-3\n" files[i] val
end
# =======================
# PLOT SETTINGS
# =======================
xticks_lab = 0:0.2:7.7

labels = ["V05612B", "V01406A", "P00748B"]
styles = [:dot, :dash, :solid]
markers = [:circle, :square, :diamond]
colors = [:orange, :deepskyblue, :purple]

# =======================
# PLOT σ
# =======================
p1 = plot()

for i in 1:3
    plot!(depth[i], σ[i],
        label=labels[i],
        linestyle=styles[i],
        lw=1.5,
        marker=markers[i],
        color=colors[i],
        xticks=xticks_lab
    )
end

plot!(xlabel="Depth [cm]",
    ylabel="σ [C cm⁻¹ s⁻¹ V⁻¹]",
    frame=:box,
    gridalpha=0.1,
    size=(600, 450))

display(p1)

# =======================
# PLOT Nd
# =======================
p2 = plot()

for i in 1:3
    plot!(depth[i], Nd[i],
        label=labels[i],
        linestyle=styles[i],
        lw=1.5,
        marker=markers[i],
        color=colors[i],
        xticks=xticks_lab
    )
end

hline!([ni],
    color=:green,
    label="nᵢ (300 K)")

plot!(xlabel="Depth [mm]",
    ylabel="Nd [cm⁻³]",
    frame=:box,
    yscale=:log10,
    ylims=(1e12, 1e15),
    gridalpha=0.1,
    size=(600, 450))

display(p2)

# =======================
# MODELLO LATERALE (σd)
# =======================
function nd_from_sigma_lat(sigma)
    mu_n = 36000
    q = 1.602e-19
    return sigma ./ (q * mu_n)
end

Nd_lat = nd_from_sigma_lat.(σd)

println("Lateral model computed.")

"""
function nd_from_sigma_lat(sigma)
    mu_n = 36000
    q = 1.602e-19
    return sigma ./ (q * mu_n)
end
Nd_lat = nd_from_sigma_lat(σd1)
"""

# σ_rad (termine dentro la radice del modello Nd)
A = 4 * mu_e * mu_h * ni^2 * e^2

depthR = [
    38 .- R[1] / 2,
    37.2 .- R[2] / 2,
    33 .- R[3] / 2
]

p3 = plot()
xticks_lab = 0:0.01:0.1
for i in 1:3
    plot!(depthR[i], σd[i],
        label=labels[i],
        linestyle=styles[i],
        lw=1.5,
        marker=markers[i],
        color=colors[i],
        #xticks=xticks_lab,
        legend=:topright
    )
end

plot!(xlabel="Depth [mm]",
    ylabel="σ\$_{rad}\$ [C cm⁻¹ s⁻¹ V⁻¹]",
    frame=:box,
    gridalpha=0.1,
    size=(600, 450))

display(p3)