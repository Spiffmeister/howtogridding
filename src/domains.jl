

function TransfiniteInterpolation(c1::Function,c2::Function,c3::Function,c4::Function,u::Real,v::Real)
    S = (1-u)*(1-v)*c1(u,v) + (1-u)*v*c2(u,v) + u*(1-v)*c3(u,v) + u*v*c4(u,v)
end








"""
"""
function GenerateDomain_Circle(R::Real,nx::Int,ny::Int)

    u = LinRange(0,1,nx)
    v = LinRange(0,1,ny)






end


function GenerateDomain_PoloidalSlice(R::Real,nx::Int,ny::Int)

    u = LinRange(0,1,nx)
    v = LinRange(0,1,ny)

end


