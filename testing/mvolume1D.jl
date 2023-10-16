
using Pkg
Pkg.activate(".")
using howtogridding


Dom1V = Grid1D([0.0,1.0],21)



D1 = Grid1D([0.0,0.5],11)
D2 = Grid1D([0.5,1.0],11)
Dom2V = GridMultiBlock(D1,D2)



D1 = Grid1D([0.0,0.35], 7)
D2 = Grid1D([0.35,0.65],7)
D3 = Grid1D([0.65,1.0], 7)
Dom3V = GridMultiBlock(D1,D2,D3)

