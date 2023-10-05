

# abstract type GridType{DIM,COORD,dtype<:Real} end
abstract type GridType{dtype<:AbstractFloat,DIM,COORD} end
abstract type LocalGridType{dtype,DIM,COORD} <: GridType{dtype,DIM,COORD} end

abstract type MetricType{MType} end


struct CartesianMetric <: MetricType{:cartesian} end
struct CurvilinearMetric <: MetricType{:curvilinear} end

Base.show(M::CartesianMetric) = print("Cartesian Metric")
Base.show(M::CurvilinearMetric) = print("Curvilinear Metric")

"""
    Grid1D
Grid data structure for 1 dimensional problems.
    
    `Grid1D{T}(𝒟::Vector,n::Integer)`
Inputs:
- Vector of domain boundaries `[x_0,x_N]`
- Number of grid points.

Returns:
- Struct for 1D grid object containing vector of grid points, ``\\Delta x`` and ``n``.
"""
struct Grid1D{TT<:Real,
        MET <:MetricType,
        GT  <:Vector{TT},
        DT  <:Union{Real,Vector{TT}}
            } <: LocalGridType{TT,1,MET}

    grid    :: GT
    Δx      :: DT
    n       :: Int64

end
function Grid1D(𝒟::Vector{TT},n::Integer) where TT
    Δx = (𝒟[2]-𝒟[1])/(n-1)
    x = collect(range(𝒟[1],𝒟[2],length=n))
    return Grid1D{TT,CartesianMetric,typeof(x),typeof(Δx)}(x,Δx,n)
end
function Grid1D(𝒟::Vector{TT}) where TT
    Δx = diff(𝒟)
    return Grid1D{TT,CurvilinearMetric,typeof(x),typeof(Δx)}(𝒟,Δx,length(𝒟))
end


"""
    Grid2D
Grid data structure for 2 dimensional problems.

    `Grid2D{T}(𝒟x::Vector,𝒟y::Vector,nx::Integer,ny::Integer)`
Inputs:
- Domain boundaries in ``x``
- Domain boundaries in ``y``
- Number of nodes in ``x``
- Number of nodes in ``y``

Returns:
- Struct for 2D grid object containing grid points in ``x`` and ``y``, ``\\Delta x`` and ``\\Delta y``, and ``n_x`` and ``n_y``.
"""
struct Grid2D{TT,
        MET<:MetricType,
        GT<:AbstractArray{TT},
        DT<:Union{Real,AbstractArray{TT}}
            } <: LocalGridType{TT,2,MET}
    gridx   :: GT
    gridy   :: GT
    Δx      :: DT
    Δy      :: DT
    nx      :: Integer
    ny      :: Integer
end
"""
    Grid2D(𝒟x::Vector,𝒟y::Vector,nx::Integer,ny::Integer)
Construct a 2D grid from the domain boundaries in ``x`` and ``y`` and the number of nodes in ``x`` and ``y``.
"""
function Grid2D(𝒟x::Vector{TT},𝒟y::Vector{TT},nx::Integer,ny::Integer) where TT
    gx = Grid1D(𝒟x,nx)
    gy = Grid1D(𝒟y,ny)
    return Grid2D{TT,CartesianMetric,typeof(gx.grid),typeof(gx.Δx)}(gx.grid, gy.grid, gx.Δx, gy.Δx, gx.n, gy.n)
end
"""
    Grid2D(𝒟x::Vector,𝒟y::Vector)
Construct a 2D grid from vectors in ``x`` and ``y``.
"""
function Grid2D(𝒟x::Vector{TT},𝒟y::Vector{TT}) where TT
    Δx = diff(𝒟x)
    Δy = diff(𝒟y)
    return Grid2D{TT,CurvilinearMetric,typeof(𝒟x),typeof(Δx)}(𝒟x,𝒟y,Δx,Δy,length(𝒟x),length(𝒟y))
end
"""
    Grid2D(cbottom::Function,cleft::Function,cright::Function,ctop::Function,nx::Integer,ny::Integer)
Construct a 2D grid from the boundary functions in ``x`` and ``y`` and the number of nodes in ``x`` and ``y``.

Curves ``c`` are parameterised by ``u`` and ``v`` where ``u`` is the coordinate in the ``x`` direction and ``v`` is the coordinate in the ``y`` direction and where ``u`` and ``v`` are in the range ``[0,1]``.
"""
function Grid2D(cbottom::Function,cleft::Function,cright::Function,ctop::Function,nx::Integer,ny::Integer)
    X,Y = meshgrid(cbottom,cleft,cright,ctop,nx,ny)

    dX = zeros(eltype(X),(nx,ny))
    dY = zeros(eltype(Y),(nx,ny))

    for i = 1:nx-1
        dX[i,:] = X[i,:] - X[i+1,:]
    end
    for j = 1:ny-1
        dY[:,j] = Y[:,j] - Y[:,j+1]
    end

    return Grid2D{eltype(X),CurvilinearMetric,typeof(X),typeof(dX)}(X, Y, dX, dY, nx, ny)
end


"""
    GridMultiBlock{TT,DIM,MET,TG,TJ,IT} <: GridType{TT,DIM,MET}
Grid data for SAT boundary problems
"""
struct GridMultiBlock{TT  <: Real,
        DIM,
        MET,
        TG,
        TJ,
        IT} <: GridType{TT,DIM,MET}
    
    Grids   :: TG
    Joint   :: TJ
    inds    :: IT
    ngrids  :: Int
end
"""
    GridMultiBlock(grids::Tuple{Vararg{Grid1D{TT,MET},N}},joints) where {N,TT,MET}
Multiblock grid for 1D grids, assumes the grids are stacked one after the other from left to right
"""
function GridMultiBlock(grids::Tuple{Vararg{Grid1D{TT,DT,MET},N}},joints) where {N,TT,DT,MET}
    inds = [sum([grids[j].n for j in 1:i]) for i in 1:length(grids)]
    return GridMultiBlock{TT, 1, MET, typeof(grids), typeof(joints),typeof(inds)}(grids,joints,inds,length(inds))
end
"""
    GridMultiBlock(grids::Tuple{Vararg{Grid2D{TT,MET},N}}) where {N,TT,MET}
Multiblock grid for 2D grids, assumes the grids are stacked in X
"""
function GridMultiBlock(grids::Tuple{Vararg{Grid2D{TT,GT,DT,MET},N}},joints) where {N,TT,GT,DT,MET}
    inds = [sum([grids[j].nx] for j in 1:i) for i in 1:length(grids)]
    return GridMultiBlock{TT,2, MET,typeof(grids),typeof(joints),typeof(inds)}(grids,joints,inds,length(inds))
end
"""
    GridMultiBlock(grids::Vector{Grid1D{TT,MET}}) where {TT,MET}
Multiblock grid for 1D grids, assumes the grids are stacked one after the other from left to right
"""
function GridMultiBlock(grids::Tuple{Vararg{Grid1D{TT,DT,MET},N}}) where {N,TT,DT,MET}
    J = [(i,i+1,Right) for i in 1:length(grids)-1]
    J = tuple(J...)
    GridMultiBlock(grids,J)
end
"""
    GridMultiBlock(grids::AbstractArray{Grid2D{TT,MET}},J::Vector{Tuple{Int,Int,NodeType}}) where {TT,MET}
Assumes blocks are stacked in ``x``
"""
function GridMultiBlock(grids::AbstractArray{Grid2D{TT,MET}}) where {TT,MET}
    J = [(i,i+1,Right) for i in 1:length(grids)-1]
    J = tuple(J...)
    GridMultiBlock(grids,J)
end



"""
    GetMinΔ
Return the miniumum grid size between ``\\Delta x`` and ``\\Delta y``
"""
function GetMinΔ end
GetMinΔ(grid::Grid1D) = grid.Δx
GetMinΔ(grid::Grid2D) = min(grid.Δx,grid.Δy)






Base.getindex(G::Grid1D,i...) = G.grid[i...]
Base.getindex(G::Grid2D,i::Integer,j::Integer) = (G.gridx[i],G.gridy[j])

function Base.getindex(G::GridMultiBlock{TT,1},i::Integer) where TT
    ii = findfirst(x->x ≥ i, G.inds)
    ii == 1 ? iii = i : iii = i - G.inds[ii-1]
    return G.Grids[ii].grid[iii]
end
function Base.getindex(G::GridMultiBlock{TT,2},i::Integer,j::Integer) where TT
    ii = findfirst(x->x ≥ i, G.inds[1,:])
    ii == 1 ? iii = i : iii = i - G.inds[1,ii-1]
    
    jj = findfirst(x->x ≥ j, G.inds)
    jj == 1 ? jjj = j : jjj = j - G.inds[jj-1]
    return G.Grids[ii,jj].grid[iii,jjj]
end




"""
    size(G::GridType)
"""
Base.size(G::Grid1D) = (G.n,)
Base.size(G::Grid2D) = (G.nx,G.ny)
function Base.size(G::GridMultiBlock{TT}) where {TT}
    # n = size.(G)
    sz = (0,0)
    for i = 1:G.ngrids
        sz = sz .+ size(G.Grids[i])
    end
    return sz
end


"""
    length(G::GridType)
"""
Base.length(G::GridType) = prod(size(G))

Base.ndims(G::GridType{TT,DIM,AT}) where {TT,DIM,AT} = DIM

Base.eachindex(G::GridType) = Base.OneTo(length(G))

eachgrid(G::GridMultiBlock) = Base.OneTo(length(G.Grids))


Base.lastindex(G::Grid1D) = G.n
Base.lastindex(G::Grid2D) = size(G)


Base.eltype(G::GridType{TT}) where TT = TT