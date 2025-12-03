cd("/home/ritaferi/Phd/SSDdev")
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using SolidStateDetectors
using Unitful
using Plots
using LegendHDF5IO
using Random

gr()


T = Float64

sim = Simulation{T}(SSD_examples[:IVCIlayer])

Random.seed!(1234) # Choose any fixed seed

mm = T(1 / 1000)

calculate_electric_potential!(sim, max_n_iterations=10, grid=Grid(sim), verbose=true, depletion_handling=true)

g = sim.electric_potential.grid

ax1, ax2, ax3 = g.axes

tick_dis = 0.02 * mm

user_additional_ticks_ax1 = sort(vcat(ax1.interval.left:tick_dis:ax1.interval.right))

user_additional_ticks_ax3 = sort(vcat(ax3.interval.left:tick_dis:ax3.interval.right))[2:end] # only support even number of ticks in z-direction 

user_ax1 = typeof(ax1)(ax1.interval, user_additional_ticks_ax1)

user_ax3 = typeof(ax3)(ax3.interval, user_additional_ticks_ax3)

user_g = typeof(g)((user_ax1, ax2, user_ax3))

calculate_electric_potential!(sim, refinement_limits=0.1, use_nthreads=8, grid=user_g, depletion_handling=true)


final_plot = plot(
    plot(sim.electric_potential, φ=0), # initial electric potential (boundary conditions)
    plot(sim.point_types), # map of different point types: fixed point / inside or outside detector volume / depleted/undepleted
    plot(sim.imp_scale), # impurity scale
    plot(sim.q_eff_imp), # charge density distribution
    layout=(1, 4), size=(3000, 1000),
    title="grid tick $max_tick_distance and V_bias $V_bias "
)
savefig(final_plot, "plots/epot_problem/claudia_plot")