using Pkg
Pkg.activate(".")  # Assicurati di essere in ~/Phd/ChargeTrappingModel_dev/SSD_code
push!(LOAD_PATH, joinpath(pwd(), "SolidStateDetectors_local_version", "src"))
using Plots
using Unitful
using SolidStateDetectors
T = Float64;

# the geometry parameters of the model for following displaydet_rin = 1u"mm"
det_z = det_r = 10u"mm"
z_draw = det_z / 2
# CHANGE THIS WITH THE SPECIFIC VALUE FOR EACH DETECTOR
pn_r = 8.957282u"mm" # this one was calculated by searching the zero impurity point (displayed in the following section)

sim = Simulation{T}(SSD_examples[:TrueCoaxial])
cfn = SSD_examples[:TrueCoaxial]
p = plot(sim.detector, xunit=u"mm", yunit=u"mm", zunit=u"mm")
display(p)


# display the impurity profile
r_list = 0u"mm":0.01u"mm":det_r
# sia pt un punto di coordinate cilindriche in cui la coordinata z è fissa, mentre la coordinata r scorre tra i differenti valori dell
# quello che fa qui è prendere la coincentrazione dei donori e degli accettori. Il modello da seguire per l'uno e per l'altrp è specificato nel file yaml.
# per il modello dei donori solitamente si usa la concentrazione con l'error function
imp_list = map(r -> let pt::CylindricalPoint{T} = CylindricalPoint(r, 0u"°", z_draw)
        SolidStateDetectors.get_impurity_density(sim.detector.semiconductor.impurity_density_model, pt) * 1e-6u"cm^-3"
    end, r_list)
#plot
p = plot(r_list, imp_list, xlabel="r / mm", ylabel="Impurity density / cm\$^{-3}\$", unitformat=:nounit, label="",
    color=:darkblue, lw=2, grid=:on, xlims=(0, 10), ylims=(-2e10, 1e11))
#  dove la concentrazione è nulla, c'è il punto di giunzione pn, che è quello che ci interessa per il nostro modello. Quindi cerchiamo quel punto e lo evidenziamo con una linea verticale.
vline!([pn_r], lw=2, ls=:dash, color=:darkred, label="PN junction boundary")
display(p)





using SolidStateDetectors: Electron, Hole
#voglio ora ottenere la mobilità degli elettroni e delle lacune in funzione della profondità, per vedere come si comportano vicino alla superficie e vicino al punto di giunzione pn. 
#Per fare questo, uso la funzione calculate_mobility che prende in input un punto e il tipo di portatore di carica (elettrone o lacuna) e restituisce la mobilità corrispondente.
#Creo una lista di punti che vanno dalla superficie (0 mm) fino al punto di giunzione pn (pn_r) e calcolo la mobilità per ciascun punto. Infine, plottiamo i risultati.
cdm = sim.detector.semiconductor.charge_drift_model
depth_list = 0u"mm":0.01u"mm":(det_r-pn_r)

# calculate_mobility calcola la mobilità di ekettorni e lacune considerando sctater con impurità ionizzate, con quelle neuutre e con i fononi
# è una funzione della profondità

mobility_list = map(depth -> let pt::CartesianPoint{T} = CartesianPoint(det_r - depth, 0, z_draw)
        (µe=SolidStateDetectors.calculate_mobility(cdm, pt, Hole) * 10000u"cm^2/(V*s)",
            µh=SolidStateDetectors.calculate_mobility(cdm, pt, Electron) * 10000u"cm^2/(V*s)")
    end, depth_list)
p = plot(depth_list, getfield.(mobility_list, :µh), label="Hole", lw=4)
plot!(depth_list, getfield.(mobility_list, :µe), label="Electron", lw=4)
plot!(xlabel="Depth to surface / mm", ylabel="Mobility / cm\$^2\$/Vs", unitformat=:nounit, legend=:topleft, xlims=(0u"mm", det_r - pn_r))
display(p)


# per calcoalre CCE dobbiamo calcolare il campo elettrico nell'inactive layer -->  per poterlo fare dobbiamo usare il codice implementato da Claudia
# dobbiamo andare a griglie mooooolto piccole nella simulazione del campo elettrico, per poter catturare bene il comportamento vicino alla superficie e vicino al punto di giunzione pn.

#La  prima cosa è calcolare il campo elettrico sulla griglia di default
calculate_electric_potential!(sim, max_n_iterations=10, grid=Grid(sim), verbose=false, depletion_handling=true)


#poi voglio andare a fare una griglia più fine, per poter catturare meglio il comportamento del campo elettrico vicino alla superficie e al punto di giunzione pn. Per fare questo, uso la funzione refine_grid che prende in input una griglia e un fattore di raffinamento e restituisce una nuova griglia più fine.
g = sim.electric_potential.grid
ax1, ax2, ax3 = g.axes

bulk_tick_dis = 0.01u"mm"
dl_tick_dis = 0.01u"mm"

user_additional_ticks_ax1 = sort(vcat(ax1.interval.left*u"m":bulk_tick_dis:pn_r, pn_r:dl_tick_dis:ax1.interval.right*u"m"))
user_ax1 = typeof(ax1)(ax1.interval, SolidStateDetectors.to_internal_units.(user_additional_ticks_ax1))
user_g = typeof(g)((user_ax1, ax2, ax3))

calculate_electric_potential!(sim, refinement_limits=0.1, grid=user_g, depletion_handling=true)
calculate_electric_field!(sim)
#calculate_weighting_potential!(sim, 1, depletion_handling = true)
#calculate_weighting_potential!(sim, 2, depletion_handling = true);

plot(
    begin
        imp = plot(sim.imp_scale, φ=0, xunit=u"mm", yunit=u"mm", title="impurity scale")
        vline!([pn_r], lw=2, ls=:dash, color=:darkred, label="PN junction boundary", legendfontsize=6)
    end,
    begin
        plot(sim.point_types, φ=0, xunit=u"mm", yunit=u"mm", title="point types")
        vline!([pn_r], lw=2, ls=:dash, color=:darkred, label="PN junction boundary", legendfontsize=6)
    end,
    begin
        plot(sim.electric_potential, xunit=u"mm", yunit=u"mm", title="electric potential")
        vline!([pn_r], lw=2, ls=:dash, color=:darkred, label="PN junction boundary", legendfontsize=6)

    end,
    begin
        plot(sim.electric_field, xunit=u"mm", yunit=u"mm", title="electric field", clims=(0, 100 * 2000))
        vline!([pn_r], lw=2, ls=:dash, color=:darkred, label="PN junction boundary", legendfontsize=6)
    end, size=(1000, 800), layout=(2, 2),
)