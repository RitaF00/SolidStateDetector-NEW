cd("/home/ritaferi/Phd/SSDdev")
using Pkg
Pkg.activate(".")

using SolidStateDetectors
using Unitful
using Plots
using Interpolations

gr()

# ho tenuto fissato la grid in modo tale da non aggiungere 
# dei punti aggiuntivi alla griglia quando cambio il max tick distance

#max_tick_distance = 5 mm
sim1 = Simulation(SSD_examples[:InvertedCoax])
sim1.detector = SolidStateDetector(sim1.detector, contact_id=2, contact_potential=500u"V")
calculate_electric_potential!(sim1,
    refinement_limits=0.2,
    verbose=false,
    depletion_handling=true,
    grid=Grid(sim1, max_tick_distance=0.5u"mm"))

# max_tick_distance = 1 mm
sim2 = Simulation(SSD_examples[:InvertedCoax])
sim2.detector = SolidStateDetector(sim2.detector, contact_id=2, contact_potential=500u"V")
calculate_electric_potential!(sim2,
    refinement_limits=0.2,
    verbose=false,
    depletion_handling=true,
    grid=Grid(sim2, max_tick_distance=0.1u"mm"))


# sono tutti SSD object

r1 = (sim1.world.intervals[1]) # radial axes (big grid)
φ1 = (sim1.world.intervals[2]) # phi axes
z1 = (sim1.world.intervals[3]) # z axes

r2 = (sim2.world.intervals[1]) # radial axes (big grid)
φ2 = (sim2.world.intervals[2]) # phi axes
z2 = (sim2.world.intervals[3]) # z axes

# phi symmetry
V1_slice = sim1.electric_potential.data[:, 1, :]; # all r, φ=0, all z
V2_slice = sim2.electric_potential.data[:, 1, :];
V1_slice = Float64.(V1_slice)
V2_slice = Float64.(V2_slice)
# la dimensione di V1 dipende dalla dimensione di r = 44 e z = 52
# la dimensione di V2 dipende dalla dimensione di r = 68 e z = 116

# r1_arr e z1_arr definiscono i punti di griglia in coordinate r e z per la matrice V1_slice.

r1_arr = Float64.(collect(range(r1.left, stop=r1.right, length=length(V1_slice[:, 1]))));
z1_arr = Float64.(collect(range(z1.left, stop=z1.right, length=length(V1_slice[1, :]))));


r2_arr = Float64.(collect(range(r2.left, stop=r2.right, length=length(V2_slice[:, 1]))));
z2_arr = Float64.(collect(range(z2.left, stop=z2.right, length=length(V2_slice[1, :]))));

itp1 = interpolate(Float64.(V1_slice), BSpline(Linear()), OnGrid());
#=
interpolate crea un interpolante B-spline lineare sulla griglia di V1_slice.

BSpline(Linear()) → interpolazione lineare tra i punti della griglia.

OnGrid() → indica che la griglia di input è regolare e corrisponde agli indici della matrice.

In pratica: itp1 è una funzione che può restituire valori interpolati di V1_slice per coordinate intermedie 
tra i punti della griglia originale.
=#

r1_range = range(r1.left, stop=r1.right, length=size(V1_slice, 1))
z1_range = range(z1.left, stop=z1.right, length=size(V1_slice, 2))

r2_range = range(r2.left, stop=r2.right, length=size(V2_slice, 1))
z2_range = range(z2.left, stop=z2.right, length=size(V2_slice, 2))


itp1_scaled = scale(itp1, r1_range, z1_range);

V1_on_fine_grid = [itp1_scaled(r2_arr[i], z2_arr[j]) for i in 1:length(r2_arr), j in 1:length(z2_arr)]

V_diff = (V2_slice - V1_on_fine_grid)  # normalizzo rispetto al potenziale di bias



diff_E_pot = heatmap(
    r2_arr,        # o r2[i].lo se non hai ancora trasformato in array
    z2_arr,        # z2[i].lo
    V_diff',
    xlabel="r [mm]",
    ylabel="z [mm]",
    colorbar_title="ΔΦ [V]",
    title="Φ small - Φ large grid",
    color=:viridis
)
plot!(sim2.detector, st=:slice, φ=0, lw=1, full_det=false, legend=false)
plot!(xlim=(minimum(r2_arr), maximum(r2_arr)), ylim=(minimum(z2_arr), maximum(z2_arr)))

savefig(diff_plot, "plots/diff_E_pot.png")