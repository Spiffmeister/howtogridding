

function plotgrid(G::GridMultiBlock)
    plt = scatter(G.Grids[1].gridx[:],G.Grids[1].gridy[:],color=:black)
    for i = 2:G.ngrids
        scatter!(G.Grids[i].gridx[:],G.Grids[i].gridy[:],color=:black)
    end
    return plt
end

