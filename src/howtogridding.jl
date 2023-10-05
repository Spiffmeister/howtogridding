module howtogridding

    using RecipesBase

    include("types.jl")
    include("grid.jl")
    include("generating.jl")

    export Grid1D, Grid2D, GridMultiBlock
    export meshgrid
    export Up, Down, Left, Right, Internal

end # module howtogridding
