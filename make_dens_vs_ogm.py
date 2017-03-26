import sys, os
from numpy import *


def Ek(kx,ky,mu):
    return -2.*1.0*(cos(kx)+cos(ky)) - mu
    
def get_dens(folder, label):
    Sig = load(folder+'Sigma.npy')

    #get mu
    with open(folder+'input'+label,'r') as f:
        for line in f:
            if "mu" in line:
                mu = float(line[line.index('#')+1:])
            if "beta" in line:
                beta = float(line[line.index('#')+1:])
                
    N_cluster = 8
    dens = 0.
    ct = 0
    Nw = 2000

    padded = zeros([8,8,Nw], dtype=complex)
    mid = shape(Sig)[2]/2
    for i in range(shape(Sig)[2]):
        padded[:,:,Nw/2-mid+i] = Sig[:,:,i,0,0]
        
    for ik1 in range(N_cluster):
        for ik2 in range(N_cluster):
            kx = -pi + 2*pi/N_cluster*ik1
            ky = -pi + 2*pi/N_cluster*ik2
            myE = Ek(kx,ky,mu)
            dens += 0.5
            for i in range(Nw):
                dens += 1./beta*1./((2.*(i-Nw/2)+1)*pi*1j/beta - myE - padded[ik1,ik2,i])
                    
    return 2.*dens/N_cluster**2
                

dirpath = sys.argv[1]

mu_map = load('mu_map.npy')
[d1,d2,d3] = shape(mu_map)

x_ogm = zeros(shape(mu_map))
d_ogm = zeros(shape(mu_map))

folders = os.listdir(dirpath)

for folder in folders:

    [i1,i2,i3] = [pos for pos,char in enumerate(folder) if char=='_']

    i = int(folder[i1+1:i2])
    j = int(folder[i2+1:i3])
    k = int(folder[i3+1:])

    label = '_%d'%i+'_%d'%j+'_%d'%k
    
    print folder
    files = os.listdir(dirpath+folder)

    d = get_dens(dirpath+folder+'/', label)

    print d

    
    d_ogm[i,j,k] = real(d)
    
       
    
print 'saving files'
save('ME_dens_ogm', d_ogm)
                      

            



