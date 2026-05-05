using CSV
using DataFrames
using Plots
using Measures

df = CSV.read(
    "diametrallygrinded.csv",
    DataFrame;
    decimal=','
)

diametro = df[!, "diameter (cm)"]
res_d = df[!, "resistance d (omh.cm)"]
lunghezza = df[!, "length (cm)"]
res_l = df[!, "resistance l (omh.cm)"]

p1 = scatter(diametro, res_d, label="")
plot!(xlabel="Diameter (cm)",
    ylabel="Resististance (Ω)",
    frame=:box,
    gridalpha=0.1,
    gridstyle=:dash)
p2 = scatter(lunghezza, res_l, label="")
plot!(xlabel="Length (cm)",
    ylabel="Resististance (Ω)",
    frame=:box,
    gridalpha=0.1,
    gridstyle=:dash)


plot(p1, p2,
    layout=(1, 2),
    size=(750, 300),
    left_margin=10Plots.mm,
    bottom_margin=10Plots.mm,
    dpi=300
)