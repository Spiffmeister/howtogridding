using Pkg
Pkg.activate(".")
using howtogridding


Dom1V = Grid2D([0.0,1.0],[0.0,1.0],21,21)



D1 = Grid2D([0.0,0.5],[0.0,1.0],11,21)
D2 = Grid2D([0.5,1.0],[0.0,1.0],11,21)

joints = (Joint(2,Right),Joint(1,Left))

Dom2V = GridMultiBlock((D1,D2),joints)