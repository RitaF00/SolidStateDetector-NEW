
using Pkg
Pkg.activate(expanduser("~/Phd/SSDdev"))
Pkg.instantiate()
using Plots
using Unitful
using SolidStateDetectors


println("miao")

x = [1, 2]
y = [1, 2]
plot(x, y)
savefig("miao.png")