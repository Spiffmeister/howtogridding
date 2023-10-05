


function coordinate(cbottom::Function,cleft::Function,cright::Function,ctop::Function,u::TT,v::TT) where TT
    S = (one(TT)-v)*cbottom(u) + v*ctop(u) + (one(TT)-u)*cleft(v) + u*cright(v) - 
        (u*v*ctop(one(TT)) + u*(one(TT)-v)*cbottom(one(TT)) + v*(one(TT)-u)*ctop(zero(TT)) + (one(TT)-u)*(one(TT)-v)*cbottom(zero(TT)))
    return S
end






"""
    meshgrid(S,TT,nx,ny)
"""
function meshgrid(TT,cbottom::Function,cleft::Function,cright::Function,ctop::Function,nx::Int,ny::Int)
    # TT = Float64

    u = LinRange(TT(0),TT(1),nx)
    v = LinRange(TT(0),TT(1),ny)

    S(u,v) = coordinate(cbottom,cleft,cright,ctop,u,v)

    
    
    X = zeros(nx,ny)
    Y = zeros(nx,ny)

    for j = 1:ny
        for i = 1:nx
            X[i,j] = S(u[i],v[j])[1]
            Y[i,j] = S(u[i],v[j])[2]
        end
    end
    return X,Y
end
"""
    meshgrid(S,nx,ny)
"""
meshgrid(cbottom::Function,cleft::Function,cright::Function,ctop::Function,nx::Int,ny::Int) = meshgrid(Float64,cbottom,cleft,cright,ctop,nx,ny)





