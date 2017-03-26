from numpy import *

N = 64
W = 8.

omegas = linspace(0.1,1.9,10)
lambdas = arange(0.02,0.46,0.04)

save('omegas',omegas)
save('lambdas',lambdas)

mu_map = zeros([len(omegas), len(lambdas), 11])

for i in range(len(omegas)):
    omega = omegas[i]
    for j in range(len(lambdas)):
        l = lambdas[j]

        g = sqrt(l * omega**2 * W / N)
        #mus = linspace(-l*W-4., -l*W+4., 11)

        mu_map[i,j,:] = linspace(-4.,0.,11)
                                                    
save('mu_map', mu_map)
print 'done saving mu map'
