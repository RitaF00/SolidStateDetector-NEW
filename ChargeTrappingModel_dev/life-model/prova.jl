using Random
using StaticArrays
using LinearAlgebra
using Plots
using ProgressMeter

"""
05/05/26

Posso generare N coppie elettrone lacuna che evolvono indipendentemente.
Rilascio puntuale di energia. --> in futuro dovrò generalizzare anche il numero di depositi di energia 

"""

# =========================
# GEOMETRIA
# =========================
const xmin = 0.0
const x_fccd = 0.1

# =========================
# FISICA
# =========================
const D = 2.89e-7   # cm^2/ns
const Δt = 1.0      # ns
const Nt = 4000
const qe = 1.602e-19

# =========================
# PARTICELLA
# =========================
mutable struct Carrier
    pos::SVector{3,Float64}
    path::Vector{SVector{3,Float64}}
    alive::Bool
    status::Symbol
end

# =========================
# INIT
# =========================
make_carrier(x0, status) = Carrier(x0, [x0], true, status)

# =========================
# DIFFUSIONE
# =========================
function diffusion_step(x)
    σ = sqrt(2 * D * Δt)
    return x + SVector(σ * randn(), σ * randn(), σ * randn())
end

# =========================
# PHYSICS
# =========================
function apply_physics!(c)
    if c.pos[1] <= xmin
        c.alive = false
    end

    if c.status == :hole && c.pos[1] >= x_fccd
        c.alive = false
    end
end

# =========================
# φ_w (FCCD model)
# =========================
phi_w(x) = exp(-x / 0.02)

# =========================
# SINGLE PAIR WAVEFORM
# =========================
function simulate_pair(x0)

    e = make_carrier(x0, :electron)
    h = make_carrier(x0, :hole)

    Q = zeros(Float64, Nt)

    for i in 2:Nt

        if e.alive
            e.pos = diffusion_step(e.pos)
            apply_physics!(e)
            push!(e.path, e.pos)
        end

        if h.alive
            h.pos = diffusion_step(h.pos)
            apply_physics!(h)
            push!(h.path, h.pos)
        end

        # Shockley–Ramo charge contribution
        Q[i] = -qe * (phi_w(e.pos[1]) - phi_w(h.pos[1]))

        if !e.alive && !h.alive
            break
        end
    end

    return Q
end

# =========================
# MULTI-PAIR EVENT
# =========================
function simulate_event(x0, Npairs)

    total_Q = zeros(Float64, Nt)

    p = Progress(Npairs, desc="Simulating event")

    for _ in 1:Npairs
        total_Q .+= simulate_pair(x0)
        next!(p)
    end

    return total_Q
end

# =========================
# RUN
# =========================
x0 = SVector(0.05, 0.0, 0.0)

E0 = 2.9          # eV per pair
E = 1000      # eV (1 keV)

Npairs = floor(Int, E / E0)

println("Generating $Npairs pairs")

waveform = simulate_event(x0, Npairs)

time = collect(0:Δt:(Nt-1)*Δt)

# =========================
# PLOT WAVEFORM
# =========================
plot(
    time,
    waveform .* 1e15,
    xlabel="time [ns]",
    ylabel="induced charge (a.u.)",
    label="Q(t) total event",
    lw=2
)