using Random
using Plots
using SpecialFunctions
using ProgressMeter
using Distributions
using JSON
using StatsBase

# -----------------------------
# Parametri fisici
# -----------------------------
Lx = 0.2
Ly = 0.2
Lz = 0.2

dx = 0.0010   # slice per diffusione Li (1 µm)
FCCD_cm = 0.1

# -----------------------------
# Annealing Li
# -----------------------------
t_ann = 18 * 60
T_ann = 623

R = 1.98
H = 11800
D0 = 2.5e-3

Ns = 10^(21.27 - 2610 / T_ann)
D_Li = D0 * exp(-H / (R * T_ann))

α = 1.6e-11

# -----------------------------
# Diffusione carica
# -----------------------------
D = 28.9
Δt = 1.0
t_max = 10000
Nt = Int(t_max / Δt)

σ = sqrt(2 * D * Δt) * 1e-4  # cm

# -----------------------------
# Raggio precipitati
# -----------------------------
r_Li = 0.002  # 20 µm 


"""
Anche in questo caso dovrei usare diverse dimensioni dei precipitati di Li.
Quello che voglio fare è introdurre una probabilità di trapping.

Da cosa dovrebbe dipendere?
- dimensione del precipitato --> più  è grande e più è probabile che vi sia trapping
- dalla profondità (?) --> questo però sembra più inserire doppiamente la dipendenza dal numero di
precipitati prodotti ad nacerta profondità dalla superficie del detector
"""

function trapping_probability()

end



