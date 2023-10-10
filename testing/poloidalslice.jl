

# # Generate the boundary curve of a D shaped poloidal slice
# function boundarycurve(R,θ)

# end



# Generate the boundary curve of a D-shaped poloidal slice
"""
    boundarycurve(R,θ,a,δ,κ)
R : major radius
θ : poloidal angle
a : triangularity parameter
δ : elongation parameter
κ : curvature parameter
"""
# Generate the boundary curve of a D-shaped poloidal slice
"""
    boundarycurve(R,θ,a,δ,κ)
R : major radius
θ : poloidal angle
a : triangularity parameter
δ : elongation parameter
κ : curvature parameter
"""
function boundarycurve(R, θ, a, δ, κ)
    # Compute the triangularity parameter
    b = δ * R

    # Compute the boundary curve
    x = R .* cos.(θ) + a .* sin.(θ)
    z = b .* sin.(θ) + κ .* sin.(2θ)
    return x, z
end


R = 1.0
theta = 0:2π/100:2π
d = 0.5

x, y = boundarycurve(R,theta,0.5,1.0,0.5)

plot(x,y,aspect_ratio=:equal)


