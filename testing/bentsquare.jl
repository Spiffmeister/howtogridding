
using Plots
using Pkg
Pkg.activate(".")
using howtogridding

cbottom(u) = [u,0]
cleft(v) = [0,v]
cright(v) = [1+2v-2v^2,v]
ctop(u) = [u,1-3u-3u^2]

# S(u,v) = (1-v)*cbottom(u) + v*ctop(u) + (1-u)*cleft(v) + u*cright(v) - ((1-u)*(1-v)*cbottom(0) + u*v*ctop(1) + (1-v)*u*cright(0) + v*(1-u)*cbottom(1))
# X,Y = meshgrid(S,x,y)


X,Y = meshgrid(cbottom,cleft,cright,ctop,11,11)




function plotgrid(X,Y)
    plt = scatter(X[1,:],Y[1,:])
    for j = 1:size(Y)[1]
        scatter!(X[:,j],Y[:,j])
    end
    for i = 1:size(X)[1]
        scatter!(X[i,:],Y[i,:])
    end
    display(plt)
end

plotgrid(X,Y)

