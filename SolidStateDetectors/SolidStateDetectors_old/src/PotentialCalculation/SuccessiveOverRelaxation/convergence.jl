function _update_till_convergence!(pcs::PotentialCalculationSetup{T,S,3},
    convergence_limit,
    via_KernelAbstractions::Bool;
    n_iterations_between_checks=500,
    depletion_handling::Val{depletion_handling_enabled}=Val{false}(),
    only2d::Val{only_2d}=Val{false}(),
    is_weighting_potential::Val{_is_weighting_potential}=Val{false}(),
    use_nthreads::Int=Base.Threads.nthreads(),
    max_n_iterations::Int=10_000, # -1 = no limit
    verbose::Bool=true
) where {T,S,depletion_handling_enabled,only_2d,_is_weighting_potential}
    backend = _ka_get_backend(pcs.potential)
    ndrange = size(pcs.potential)[1:3] .- 2
    kernel = get_sor_kernel(S, backend, Val(via_KernelAbstractions))

    c_limit = _is_weighting_potential ? convergence_limit :
              abs(convergence_limit * (iszero(pcs.bias_voltage) ? maximum(abs.(pcs.potential)) : pcs.bias_voltage))

    # one() : resituisce l'unità moltiplicatica. se c_lim è un INT da 1, se un FLOAT64 da 1.0     
    c = (one(c_limit) + c_limit) * 10
    n_performed_iterations = 0
    tmp_potential = similar(pcs.potential, ndrange)
    inner_ranges = broadcast(i -> 2:size(tmp_potential, i)+1, (1, 2, 3))
    is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")
    cs = fill(c, 4)

    prog = nothing
    #Convergence:  (thresh = <c_limit>, value = <c>)
    if verbose
        prog = ProgressThresh(c_limit; dt=0.1, desc="Convergence: ", output=stderr, enabled=!is_logging(stderr))
    end

    stop_reason = nothing  # will hold :converged, :stable_variation, :max_iterations

    while c > c_limit
        for _ in 1:(n_iterations_between_checks-1)
            # !UPDATE è la funzione in cui si applica il SOR 
            update!(pcs, kernel, ndrange; use_nthreads, depletion_handling, is_weighting_potential, only2d)
            n_performed_iterations += 1
        end

        # one more update + misura della differenza
        # copio la versione corrente del potenziale pe rpoter fare un'ulteriore iterazione e poi la dofferenza e fare il check della convergenza.
        tmp_potential[:, :, :] .= view(pcs.potential, inner_ranges..., 1)
        update!(pcs, kernel, ndrange; use_nthreads, depletion_handling, is_weighting_potential, only2d)
        tmp_potential[:, :, :] .-= view(pcs.potential, inner_ranges..., 1)
        n_performed_iterations += 1

        c = maximum(abs.(tmp_potential))
        if verbose && prog !== nothing
            ProgressMeter.update!(prog, c)
        end

        cs = circshift(cs, -1)
        cs[end] = c
        cs_μ = mean(cs)
        cs_σ = std(cs, mean=cs_μ)

        # 1) stop per variazione stabile
        if cs_σ < c_limit
            stop_reason = :stable_variation
            println("🚫 STOP: Convergence limit not reached but value of c is stable (cs_σ < c_limit).")
            break
        end

        # 2) stop per convergenza (lo rilevo immediatamente)
        if c <= c_limit
            stop_reason = :converged
            println("🟢 STOP: convergence reached (c <= c_limit).")
            break
        end

        # 3) controllo max iterazioni
        if max_n_iterations != -1 && n_performed_iterations >= max_n_iterations
            stop_reason = :max_iterations
            println("⛔ STOP: reached maximum iteration number ($n_performed_iterations >= $max_n_iterations).")
            break
        end
    end

    if verbose && prog !== nothing
        ProgressMeter.finish!(prog)
    end

    # Se il ciclo è uscito naturalmente perché while c > c_limit diventato false
    # (cioè non è stato fatto break ma la condizione non è più vera), registriamo convergenza.
    if stop_reason === nothing
        if c <= c_limit
            stop_reason = :converged
            println("🟢 STOP: convergence reached (outside c <= c_limit).")
        elseif max_n_iterations != -1 && n_performed_iterations >= max_n_iterations
            stop_reason = :max_iterations
            println("⛔ STOP: reached maximum iteration number (post-check).")
        else
            # caso improbabile, fallback
            stop_reason = :unknown
            println("ℹ️ STOP: Killed cycle (reason unknown). c = $c")
        end
    end

    println("🔄 Number of iterations: $n_performed_iterations; reason: $stop_reason; final c = $c")
    return c
end
