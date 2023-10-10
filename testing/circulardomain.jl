#=

=#
using Plots
using Pkg

Pkg.activate(".")
using howtogridding



R = 1.0

D = GenerateDomain_Circle(R,5,5)

plt = plotgrid(D)
plt



#=
# First build the full circle for plotting
x = [R*cos(θ) for θ in 0:0.01:2π]
y = [R*sin(θ) for θ in 0:0.01:2π]




# Now the square
# square = Grid2D([-0.25,0.25],[-0.25,0.25],5.0,5.0)


# each of the square edges
squpper(u)  = [u/2 - 0.25,  0.25]
sqleft(v)   = [-0.25,       v/2 - 0.25]
sqright(v)  = [0.25,        v/2 - 0.25]
sqbottom(u)    = [u/2 - 0.25,  -0.25]


# Top annulus
AUbottom(u) = [u/2 - 0.25, 0.25]
AUleft(v)   = v*[cos(3π/4) + 0.25, sin(3π/4) - 0.25] + [-0.25, 0.25]
AUright(v)  = v*[cos(π/4) - 0.25, sin(π/4) - 0.25] + [0.25, 0.25]
AUtop(u)    = [cos(u*(π/4 - 3π/4) + 3π/4), sin(u*(π/4 - 3π/4) + 3π/4)]

# Left annulus
ALbottom(u) = u*[-cos(3π/4) + 0.25, -sin(3π/4) + 0.25] + [-0.25, 0.25]
ALleft(v)   = [cos(v*(3π/4 - 5π/4) + 5π/4), sin(v*(3π/4 - 5π/4) + 5π/4)]
ALright(v)  = [-0.25, v/2 - 0.25]
ALtop(u)    = u*[0.25 + cos(3π/4), 0.25 + sin(π/4)] + [-cos(3π/4), -sin(π/4)]

# Bottom annulus
ABbottom(u) = [cos(u*(7π/4 - 5π/4) + 5π/4), sin(u*(7π/4 - 5π/4) + 5π/4)]
ABleft(v)   = v*[-0.25 - cos(5π/4), -0.25 - sin(5π/4)] + [cos(5π/4), sin(5π/4)]
ABright(v)  = v*[0.25 - cos(7π/4), -0.25 - sin(7π/4)] + [cos(7π/4), sin(7π/4)]
ABtop(u)    = [u/2 - 0.25, -0.25]

# Right annulus
ARbottom(u) = u*[cos(7π/4) - 0.25, sin(7π/4) + 0.25] + [0.25, -0.25]
ARleft(v)   = [0.25, v/2 - 0.25]
ARright(v)  = [cos(v*(9π/4 - 7π/4) + 7π/4), sin(v*(9π/4 - 7π/4) + 7π/4)]
ARtop(u)    = u*[cos(π/4) - 0.25, sin(π/4) - 0.25] + [0.25, 0.25]



nx = ny = 10

sq = Grid2D(sqbottom,sqleft,sqright,squpper,nx,ny)
AU = Grid2D(AUbottom,AUleft,AUright,AUtop,nx,ny)
AL = Grid2D(ALbottom,ALleft,ALright,ALtop,nx,ny)
AB = Grid2D(ABbottom,ABleft,ABright,ABtop,nx,ny)
AR = Grid2D(ARbottom,ARleft,ARright,ARtop,nx,ny)


glayout = ([(2,Up),(3,Left),(4,Down),(5,Right)],
            [(1,Down),(3,Left),(5,Right)],
            [(1,Right),(2,Up),(4,Down)],
            [(1,Up),(3,Right),(5,Right)],
            [(1,Left),(2,Left),(4,Down)])


G = GridMultiBlock((sq,AU,AL,AB,AR),glayout)



scatter(sq.gridx[:],sq.gridy[:],label="centre",color=:black)
scatter!(AU.gridx[:],AU.gridy[:],label="upper",color=:red)
scatter!(AL.gridx[:],AL.gridy[:],label="left",color=:blue)
scatter!(AB.gridx[:],AB.gridy[:],label="right",color=:green)
scatter!(AR.gridx[:],AR.gridy[:],label="down",color=:yellow)



# AUright(v) = v/sqtocirc .* [cos(3π/4),sin(3π/4)] + [0.25,0.25]



"""
    circlesection(x::Real,y::Real,r::Real,theta::Real)
Parameterise the section of a circles circumference with radius ``r`` and angle ``theta`` centered at ``(x,y)``.
"""
function circlesection(x::Real,y::Real,r::Real,theta::Real)
    [x+r*cos(theta),y+r*sin(theta)]
end
=#









