module howtogridding

    using Plots

    include("types.jl")
    include("basis.jl")
    include("grid.jl")
    include("generating.jl")
    include("domains.jl")
    include("plotgrid.jl")

    export Grid1D, Grid2D, GridMultiBlock
    export meshgrid
    export Up, Down, Left, Right, Internal

    export GenerateDomain_Circle
    export Torus

    export plotgrid


end # module howtogridding
