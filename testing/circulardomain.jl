#=

=#
using Plots
using Pkg

Pkg.activate(".")
using howtogridding



R = 1.0

D = GenerateDomain_Circle(R,5,5)

plt = plotgrid(D)
plt











