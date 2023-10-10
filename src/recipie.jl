# RecipesBase.@recipe function Grid(::Type{GT},G::GT)
    
#     X,Y = zeros(G.nx,G.ny),zeros(G.nx,G.ny)

# end


function plotgrid(G::GridMultiBlock)
    scatter(G.Grids[1].gridx,G.Grids[1].gridy,color=:black)
    for i = 2:G.ngrids
        scatter!(G.Grids[i].gridx,G.Grids[i].gridy,color=:black)
    end
end


