
using Pkg
using Unitful
using Plots
using HDF5
using ProgressMeter
using Random
using SpecialFunctions
using Distributions



function multiple_charges_trapping(x_charges, N::Int=100)
    # se x_charges è un singolo numero, trasformalo in array
    if isa(x_charges, Number)
        x_charges = fill(x_charges, N)
    elseif length(x_charges) != N
        error("La lunghezza di x_charges deve essere uguale a N")
    end

    # --------------------------
    # Parametri Li layer
    # --------------------------
    Lx = 0.12      # profondità in cm
    Ly = 0.2       # estensione y in cm
    dx = 0.001     # passo per generare slice
    α = 1e-8       # thinning factor

    t_ann = 18 * 60 # 18 minutes in seconds
    T_ann = 623
    Na = 3 * 10^9
    R = 1.98

    H = if 473 <= T_ann <= 873
        11800
    else
        10700
    end

    D0 = if 473 <= T_ann <= 873
        2.5 * 10^(-3)
    else
        1.3 * 10^(-3)
    end

    Ns = 10^(21.27 - 2610 / T_ann)


    FCCD_cm = 0.09 # profondità FCCD in cm
    PN = 0.11

    # Parametri diffusione Li
    H = 11800
    D0 = 2.5e-3
    Ns = 10^(21.27 - 2610 / T_ann)

    NGe = 4.42e22
    D_Li = D0 * exp(-H / (R * T_ann))
    #Nd(x) = Ns .* erfc.(x ./ (2 * sqrt(D_Li * t_ann))) ./ NGe
    Nd(x) = Ns .* erfc.(x ./ (2 * sqrt(D_Li * t_ann)))
    #=
    function Nd(x)
        c0 = 0.005       # 0.5 % come frazione
        # L_max = 0.1      # profondità massima in cm
        L_max = 0.06 # dalla tesi
        return  c0 * max(0, 1 - x / L_max)
    end
    =#
    # --------------------------
    # Generazione punti Li 2D
    # --------------------------
    x_slices = 0:dx:Lx
    x_positions = Float64[]
    y_positions = Float64[]
    for xi in x_slices
        Nd_val = Nd(xi)
        dA = dx * Ly * Ly # devo considerare la 3 dimensione
        λ = α * Nd_val * dA
        N_i = rand(Poisson(λ))
        for _ in 1:N_i
            push!(x_positions, xi + rand() * dx)
            push!(y_positions, rand() * Ly)
        end
    end

    # --------------------------
    # Parametri diffusione singola carica
    # --------------------------
    D = 28.9           # diffusione [µm^2/ns]  --> da cambiare --> mettere dipedenda da μ
    Δt = 1.0           # passo temporale [ns]
    t_max = 10000      # ns
    Nt = Int(t_max / Δt)
    σ = sqrt(2 * D * Δt)   # µm^2/ns -> µm

    # --------------------------
    # Generazione posizioni iniziali N cariche
    # --------------------------
    y_charges = rand(N) .* Ly   # posizioni y casuali tra 0 e Ly

    # --------------------------
    # Tracking raccolta
    # --------------------------
    collected_count = 0

    for n in 1:N
        x_charge = x_charges[n]
        y_charge = y_charges[n]
        trapped = false

        # se nasce oltre Lx, raccolta immediata
        if x_charge >= Lx
            collected_count += 1
            continue
        end

        # simulazione diffusione
        for i in 1:Nt
            x_charge += σ * randn() * 1e-4    # µm -> cm
            y_charge += σ * randn() * 1e-4    # µm -> cm
            x_charge = clamp(x_charge, 0, Lx)
            y_charge = clamp(y_charge, 0, Ly)

            r_Li = 1e-4 # 1 µm
            # trapping Li
            for (xi, yi) in zip(x_positions, y_positions)
                if abs(x_charge - xi) < r_Li && abs(y_charge - yi) < r_Li
                    trapped = true
                    break
                end
            end

            if trapped
                break
            end

            # raccolta FCCD
            if x_charge >= FCCD_cm
                collected_count += 1
                break
            end
        end
    end

    #println("Numero di cariche raccolte alla FCCD: $collected_count su $N")
    return collected_count
end


#---------- starting code ----------------
println("miao")
x_pos = 0:0.01:0.15
N_charges = 100
N_repeat = 1  # quante simulazioni per profondità

# matrice per salvare tutti i valori
N_matrix = zeros(Int, length(x_pos), N_repeat)

p = Progress(length(x_pos) * N_repeat, desc="Simulazione CCE")

for (i, x) in enumerate(x_pos)
    for j in 1:N_repeat
        N_matrix[i, j] = multiple_charges_trapping(x, N_charges)
        next!(p)
    end
end

CCE = N_matrix ./ N_charges

CCE_mean = vec(mean(CCE, dims=2))
CCE_std = vec(std(CCE, dims=2))


blue = palette(:auto)[1]

p = plot(
    xlabel="depth (mm)",
    ylabel="CCE",
    frame=:box,
    legend=false
)

# singole simulazioni
for j in 1:N_repeat
    plot!(
        p,
        x_pos,
        CCE[:, j],
        marker=:o,
        seriestype=:scatter,
        color=blue,
        ms=2
    )
end

# media arancione con errore
plot!(
    p,
    x_pos,
    CCE_mean,
    color=:orange,)

vline!(p, [0.11], ls=:dash, color=:green, label="P-N junction")