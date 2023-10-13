
using Pkg
Pkg.activate(".")
using howtogridding
using LinearAlgebra
using Test

xy = collect(LinRange(0.0,1.0,5))

X,Y = meshgrid(xy,xy)

D = Grid2D(X,Y)



@test repeat(collect(0.0:0.25:1.0),outer=(1,5)) == D.gridx
@test repeat(collect(0.0:0.25:1.0),outer=(1,5))' == D.gridy








cbottom(u) = [u,0.0]
cleft(v) = [0.0,v]
cright(v) = [1.0,v]
ctop(u) = [u,1.0]

D = Grid2D(cbottom,cleft,cright,ctop,5,5)



@test repeat(collect(0.0:0.25:1.0),outer=(1,5)) == D.gridx
@test repeat(collect(0.0:0.25:1.0),outer=(1,5))' == D.gridy



"""
    Testing 
"""
cbottom(u) = [u,-u/2]
cleft(v) = [0.0,v]
cright(v) = [1.0,v-1/2]
ctop(u) = [u,1.0-u/2]

D = Grid2D(cbottom,cleft,cright,ctop,5,5)


