
using Pkg
Pkg.activate(".")
using howtogridding
using LinearAlgebra

xy = collect(LinRange(0.0,1.0,5))

X,Y = meshgrid(xy,xy)

D = Grid2D(X,Y)



