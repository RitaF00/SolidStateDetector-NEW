using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Plots
using Unitful
using SolidStateDetectors

using SolidStateDetectors: Electron, Hole
T = Float64;
det_rin = 1u"mm"
det_z = det_r = 10u"mm"
z_draw = det_z / 2
pn_r = 8.957282u"mm" # this one was calculated by searching the zero impurity point (displayed in the following section)

sim = Simulation{T}(SSD_examples[:TrueCoaxial])
cfn = SSD_examples[:TrueCoaxial]
#print(open(f -> read(f, String), cfn))
p = plot(sim.detector, xunit=u"mm", yunit=u"mm", zunit=u"mm")
display(p)
r_list = 0u"mm":0.01u"mm":det_r
imp_list = map(r -> let pt::CylindricalPoint{T} = CylindricalPoint(r, 0u"°", z_draw)
        SolidStateDetectors.get_impurity_density(sim.detector.semiconductor.impurity_density_model, pt) * 1e-6u"cm^-3"
    end, r_list)
p = plot(r_list, imp_list, xlabel="r / mm", ylabel="Impurity density / cm\$^{-3}\$", unitformat=:nounit, label="",
    color=:darkblue, lw=2, grid=:on, xlims=(0, 10), ylims=(-2e10, 1e11))
vline!([pn_r], lw=2, ls=:dash, color=:darkred, label="PN junction boundary")
display(p)
using SolidStateDetectors: Electron, Hole
cdm = sim.detector.semiconductor.charge_drift_model
depth_list = 0u"mm":0.01u"mm":(det_r-pn_r)
mobility_list = map(depth -> let pt::CartesianPoint{T} = CartesianPoint(det_r - depth, 0, z_draw)
        (µe=SolidStateDetectors.calculate_mobility(cdm, pt, Hole) * 10000u"cm^2/(V*s)",
            µh=SolidStateDetectors.calculate_mobility(cdm, pt, Electron) * 10000u"cm^2/(V*s)")
    end, depth_list)

p = plot(depth_list, getfield.(mobility_list, :µh), label="Hole", lw=4)
plot!(depth_list, getfield.(mobility_list, :µe), label="Electron", lw=4)
plot!(xlabel="Depth to surface / mm", ylabel="Mobility / cm\$^2\$/Vs", unitformat=:nounit, legend=:topleft, xlims=(0u"mm", det_r - pn_r))
display(p)