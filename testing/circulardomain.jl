
using Plots




theta(x,y) = atan(y/x)


"""
    circlesection(x::Real,y::Real,r::Real,theta::Real)
Parameterise the section of a circles circumference with radius ``r`` and angle ``theta`` centered at ``(x,y)``.
"""
function circlesection(x::Real,y::Real,r::Real,theta::Real)
    [x+r*cos(theta),y+r*sin(theta)]
end




"""
    Annulus(verticies::Vector{Vector{Real}})
verticies is a vector of vectors of the form [x,y]
"""
function Annulus(verticies::Vector{Vector{Real}})

    # Left edge
    verticies[1][1]*x + verticies[2][1]*x



    vertex[1.0]

end






inner = Grid2D([-1.0,1.0],[-1.0,1.0],5,5)



