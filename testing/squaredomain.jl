
using Plots


c1(u) = [(2u-1), -1]
c2(u) = [(2u-1), 1]
c3(v) = [1, (2v-1)]
c4(v) = [-1,(2v-1)]


S(u,v) = (1-v)*c1(u) + v*c3(u) + (1-u)*c2(v) + u*c4(v) - ((1-u)*(1-v)*c1(0) + u*v*c3(1) + (1-v)*u*c4(0) + v*(1-u)*c1(1))


x = LinRange(0.0,1.0,11)
y = LinRange(0.0,1.0,11)




X = zeros(length(x),length(y))
Y = zeros(length(x),length(y))

for j = 1:length(y)
    for i = 1:length(x)
        X[i,j] = S(x[i],y[j])[1]
        Y[i,j] = S(x[i],y[j])[2]
    end
end



scatter(X[:,1],Y[1,:])
for j = 1:length(y)
    scatter!(X[:,j],Y[j,:])
end

