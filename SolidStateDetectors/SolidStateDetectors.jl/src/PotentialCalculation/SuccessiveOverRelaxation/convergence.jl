# convergence.jl

# ==============================
# Pacchetti necessari
# ==============================
using JSON
using Base.Threads
using ProgressMeter
# altri pacchetti necessari dal tuo progetto
# using GPUArrays
# using SolidStateDetectors.Utils
# etc.




# =============0 funzione aggiuntiva per cambiare nome al dizionario JSON ==============
function next_available_filename(base::String)
    if !isfile(base)
        return base
    end

    name, ext = splitext(base)
    i = 1
    while true
        candidate = "$(name)_$(i)$(ext)"
        if !isfile(candidate)
            return candidate
        end
        i += 1
    end
end

function _update_till_convergence!(pcs::PotentialCalculationSetup{T,S,3},
    convergence_limit,
    via_KernelAbstractions::Bool;
    n_iterations_between_checks=500,
    depletion_handling::Val{depletion_handling_enabled}=Val{false}(),
    only2d::Val{only_2d}=Val{false}(),
    is_weighting_potential::Val{_is_weighting_potential}=Val{false}(),
    use_nthreads::Int=Base.Threads.nthreads(),
    max_n_iterations::Int=10_000,
    verbose::Bool=true
) where {T,S,depletion_handling_enabled,only_2d,_is_weighting_potential}

    backend = _ka_get_backend(pcs.potential)
    ndrange = size(pcs.potential)[1:3] .- 2
    kernel = get_sor_kernel(S, backend, Val(via_KernelAbstractions))

    c_limit = _is_weighting_potential ? convergence_limit :
              abs(convergence_limit * (iszero(pcs.bias_voltage) ?
                                       maximum(abs.(pcs.potential)) :
                                       pcs.bias_voltage))

    c = (one(c_limit) + c_limit) * 10
    n_performed_iterations = 0

    tmp_potential = similar(pcs.potential, ndrange)
    inner_ranges = broadcast(i -> 2:size(tmp_potential, i)+1, (1, 2, 3))
    cs = fill(c, 4)

    if verbose
        prog = ProgressThresh(c_limit; dt=0.1, desc="Convergence: ", output=stderr)
    end

    c_single = Float64[]
    c_array = Float64[]
    n_array = Int[]

    # 🔹 MOTIVO DI STOP
    stop_reason = "unknown"

    while c > c_limit
        for _ in 1:n_iterations_between_checks-1
            old_potential = copy(Array(pcs.potential))

            update!(pcs, kernel, ndrange;
                use_nthreads,
                depletion_handling,
                is_weighting_potential,
                only2d)

            new_potential = Array(pcs.potential)
            c_local = maximum(abs.(new_potential .- old_potential))
            push!(c_single, c_local)

            n_performed_iterations += 1
        end

        tmp_potential[:, :, :] .= view(pcs.potential, inner_ranges..., 1)
        update!(pcs, kernel, ndrange;
            use_nthreads,
            depletion_handling,
            is_weighting_potential,
            only2d)
        tmp_potential[:, :, :] .-= view(pcs.potential, inner_ranges..., 1)

        n_performed_iterations += 1
        c = maximum(abs.(tmp_potential))

        push!(c_array, c)
        push!(n_array, n_performed_iterations ÷ n_iterations_between_checks)

        if verbose
            ProgressMeter.update!(prog, c)
        end

        cs = circshift(cs, -1)
        cs[end] = c

        # ✅ STOP: convergenza stabile
        if std(cs) < c_limit
            stop_reason = "std(cs) < c_limit (stable convergence)"
            break
        end

        # ✅ STOP: massimo numero di iterazioni
        if max_n_iterations != -1 && n_performed_iterations >= max_n_iterations
            stop_reason = "max_n_iterations reached"
            break
        end
    end

    # ✅ STOP: convergenza diretta
    if c <= c_limit
        stop_reason = "c <= c_limit (direct convergence)"
    end

    if verbose
        ProgressMeter.finish!(prog)
    end

    println("🔄 Number of iterations: $n_performed_iterations; final c = $c")
    println("🛑 Loop stopped because: $stop_reason")

    # qui posso salvare l'andamento dell'errore di convergenza su file JSON
    #===
    # 🔹 Salvataggio NON distruttivo 🔹
    if _is_weighting_potential
        filename = next_available_filename("c_single.json")
        open(filename, "w") do io
            JSON.print(io, Dict(
                "c_single" => c_single,
                "stop_reason" => stop_reason,
                "final_c" => c,
                "iterations" => n_performed_iterations
            ))
        end
        println("📁 Saved c_single to: $filename")
    end

    ===#

    return c
end