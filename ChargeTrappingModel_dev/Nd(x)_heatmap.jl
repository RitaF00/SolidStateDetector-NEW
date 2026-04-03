using Plots
using SpecialFunctions

gr()    # interattivo (puoi ruotare il grafico)

# -----------------------------
# Costanti
# -----------------------------
t_ann = 18 * 60
R = 1.98
H = 11800
D0 = 2.5e-3

# -----------------------------
# Range
# -----------------------------
T_range = range(500, 900, length=120)   # Temperatura [K]
x_range = range(0, 1e-4, length=120)    # Depth [cm]

# -----------------------------
# Funzioni
# -----------------------------
Ns(T) = 10^(21.27 - 2610 / T)
D(T) = D0 * exp(-H / (R * T))

Nd(T, x) = Ns(T) * erfc(x / (2 * sqrt(D(T) * t_ann)))

# -----------------------------
# Griglia dati
# -----------------------------
Z = [Nd(T, x) for x in x_range, T in T_range]

# -----------------------------
# 3D SURFACE (normale)
# -----------------------------
p1 = surface(
    T_range,
    x_range,
    Z,
    xlabel="Temperatura (K)",
    ylabel="Depth (cm)",
    zlabel="Nd",
    title="Nd(T, x)",
    color=:viridis,
    linewidth=0
)

# -----------------------------
# 3D SURFACE LOG (consigliata)
# -----------------------------
p2 = surface(
    T_range,
    x_range,
    log10.(Z),
    xlabel="Temperatura (K)",
    ylabel="Depth (cm)",
    zlabel="log10(Nd)",
    title="Nd(T, x) - scala log",
    color=:plasma,
    linewidth=0
)

# -----------------------------
# MOSTRA (una alla volta o entrambe)
# -----------------------------
display(p1)
display(p2)

# -----------------------------
# SALVA IMMAGINI
# -----------------------------
savefig(p1, "Nd_surface.png")
savefig(p2, "Nd_surface_log.png")

# -----------------------------
# SALVA VERSIONE INTERATTIVA HTML
# -----------------------------
savefig(p1, "Nd_surface.html")
savefig(p2, "Nd_surface_log.html")