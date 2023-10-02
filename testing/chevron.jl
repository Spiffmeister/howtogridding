
using Pkg
Pkg.activate(".")
using howtogridding


function cbottom(u)
    x = u
    u ≤ 0.5 ? y = -u : y = u-1
    return [x,y]
end

cleft(v) = [0.0,v]

function ctop(u)
    x = u
    u ≤ 0.5 ? y = 1-u : y = u
    return [x,y]
end

cright(v) = [1.0,v]


X,Y = meshgrid(cbottom,cleft,cright,ctop,10,10)


using Plots

function plotgrid(X,Y)
    plt = plot(legend=false)
    for j = 1:size(Y)[1]
        plot!(X[:,j],Y[:,j],color=:black)
    end
    for i = 1:size(X)[1]
        plot!(X[i,:],Y[i,:],color=:black)
    end
    display(plt)
end

plotgrid(X,Y)





