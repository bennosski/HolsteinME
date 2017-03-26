import sys, os
from scipy.optimize import curve_fit
from numpy import *

def get_Nws_xs(folder):
    files = os.listdir(folder)
    xs = []
    Nws = []
    for myfile in files:
        if 'xsc_' in myfile:
            index = myfile.index("w")
            Nws.append(float(myfile[index+1:]))
            xs.append(loadtxt(folder+myfile))
    Nws = asarray(Nws)
    xs = abs(asarray(xs))
    
    return Nws,xs

def lin(x,a,b):
        return a*x + b

def get_xsc(folder):
    Nws,xs = get_Nws_xs(folder)

    if os.path.exists(folder+'xsc_Nw40'):
        dens = float(loadtxt(folder+'dens'))
        p,c = curve_fit(lin, 1./Nws, xs, [-0.01,0.196])
        return dens,p[1]
    else:
        print 'no data in ',folder
        return 0,0
    
dirpath = sys.argv[1]
#dirpath = 'outputfiles_2_24_17_moreBeta/'

mu_map = load('mu_map.npy')
[d1,d2,d3] = shape(mu_map)

x_ogm = zeros(shape(mu_map))
d_ogm = zeros(shape(mu_map))

folders = os.listdir(dirpath)


ct = 0
for folder in folders:

    [i1,i2,i3] = [pos for pos,char in enumerate(folder) if char=='_']

    i = folder[i1+1:i2]
    j = folder[i2+1:i3]
    k = folder[i3+1:]
            
    #print folder
    files = os.listdir(dirpath+folder)

    print folder,len(files)

    if len(files)<8:
        ct += 1

print 'remaining to be run ',ct

            



