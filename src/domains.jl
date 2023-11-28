

function TransfiniteInterpolation(c1::Function,c2::Function,c3::Function,c4::Function,u::Real,v::Real)
    S = (1-u)*(1-v)*c1(u,v) + (1-u)*v*c2(u,v) + u*(1-v)*c3(u,v) + u*v*c4(u,v)
end







"""
"""
function GenerateDomain_Circle(R::Real,nx::Int,ny::Int)

    u = LinRange(0,1,nx)
    v = LinRange(0,1,ny)


    sq_corner = R*0.25

    # each of the square edges
    # sqbottom(u) = [u*(sq_corner + sq_corner) - sq_corner,  -sq_corner]
    # sqleft(v)   = [-sq_corner,       v*(sq_corner + sq_corner) - sq_corner]
    # sqright(v)  = [sq_corner,        v*(sq_corner + sq_corner) - sq_corner]
    # squpper(u)  = [u*(sq_corner + sq_corner) - sq_corner,  sq_corner]


    squ = repeat(collect(range(-sq_corner,sq_corner,length=nx)),1,ny)
    sqv = repeat(collect(range(-sq_corner,sq_corner,length=ny))',nx,1)
    sq = Grid2D(squ,sqv)


    # Top annulus
    AUbottom(u) = [u*(sq_corner + sq_corner) - sq_corner, sq_corner]
    AUleft(v)   = v*[R*cos(3π/4) + sq_corner, R*sin(3π/4) - sq_corner] + [-sq_corner, sq_corner]
    AUright(v)  = v*[R*cos(π/4) - sq_corner, R*sin(π/4) - sq_corner] + [sq_corner, sq_corner]
    AUtop(u)    = [R*cos(u*(π/4 - 3π/4) + 3π/4), R*sin(u*(π/4 - 3π/4) + 3π/4)]

    # Left annulus
    ALbottom(u) = u*[-R*cos(3π/4) + sq_corner, -R*sin(3π/4) + sq_corner] + [-sq_corner, sq_corner]
    ALleft(v)   = [R*cos(v*(3π/4 - 5π/4) + 5π/4), R*sin(v*(3π/4 - 5π/4) + 5π/4)]
    ALright(v)  = [-sq_corner, v*(sq_corner + sq_corner) - sq_corner]
    ALtop(u)    = u*[sq_corner + R*cos(3π/4), sq_corner + R*sin(π/4)] + [-R*cos(3π/4), -R*sin(π/4)]

    # Bottom annulus
    ABbottom(u) = [R*cos(u*(7π/4 - 5π/4) + 5π/4), R*sin(u*(7π/4 - 5π/4) + 5π/4)]
    ABleft(v)   = v*[-sq_corner - R*cos(5π/4), -sq_corner - R*sin(5π/4)] + [R*cos(5π/4), R*sin(5π/4)]
    ABright(v)  = v*[sq_corner - R*cos(7π/4), -sq_corner - R*sin(7π/4)] + [R*cos(7π/4), R*sin(7π/4)]
    ABtop(u)    = [u*(sq_corner + sq_corner) - sq_corner, -sq_corner]

    # Right annulus
    ARbottom(u) = u*[R*cos(7π/4) - sq_corner, R*sin(7π/4) + sq_corner] + [sq_corner, -sq_corner]
    ARleft(v)   = [sq_corner, v*(sq_corner + sq_corner) - sq_corner]
    ARright(v)  = [R*cos(v*(9π/4 - 7π/4) + 7π/4), R*sin(v*(9π/4 - 7π/4) + 7π/4)]
    ARtop(u)    = u*[R*cos(π/4) - sq_corner, R*sin(π/4) - sq_corner] + [sq_corner, sq_corner]


    AU = Grid2D(AUbottom,AUleft,AUright,AUtop,nx,ny)
    AL = Grid2D(ALbottom,ALleft,ALright,ALtop,nx,ny)
    AB = Grid2D(ABbottom,ABleft,ABright,ABtop,nx,ny)
    AR = Grid2D(ARbottom,ARleft,ARright,ARtop,nx,ny)


    glayout = ([(2,Up),(3,Left),(4,Down),(5,Right)],
                [(1,Down),(3,Left),(5,Right)],
                [(1,Right),(2,Up),(4,Down)],
                [(1,Up),(3,Right),(5,Right)],
                [(1,Left),(2,Left),(4,Down)])


    cgrid = GridMultiBlock((sq,AU,AL,AB,AR),glayout)
    
    return cgrid

end


