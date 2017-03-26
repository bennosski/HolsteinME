from numpy import *
import subprocess, os, time, sys

##first get a sense of where mu is with sdev
#apparently mu should be around -lambda*W = -N*g^2/omega^2  where N is the total N

N = 64  # Nx**2
W = 8.

omegas = load('omegas.npy')
lambdas = load('lambdas.npy')
mu_map = load('mu_map.npy')
[d1,d2,d3] = shape(mu_map)


#dirpath = 'inputfiles_2_24_17_moreBeta/'
#dirpath = 'inputfiles_2_25_17_w4p0_l0p2/'
dirpath = sys.argv[1]


def bash_command(cmd):
    subprocess.Popen(['/bin/bash', '-c', cmd])

istart = 0
jstart = 0
kstart = 0

for i in range(len(omegas)):
    omega = omegas[i]
    for j in range(len(lambdas)):
        lamb = lambdas[j]

        g = sqrt(lamb * omega**2 * W / N) 

        for k in range(shape(mu_map)[2]):
            
            if k + j*d3 + i*d3*d2 < kstart + jstart*d3 + istart*d3*d2:
                continue
            
            mu = mu_map[i,j,k]

            label = '_%d'%i+'_%d'%j+'_%d'%k
            
            input_file_name    = dirpath+'input'+label
            
            cmd = " cp input temp1"+label
            cmd += """; sed "/g_dqmc/c\g_dqmc  # %1.15f"""%g+'" temp1'+label+' > temp2'+label
            cmd += '; cp temp2'+label+' '+'temp1'+label
            cmd +=  """; sed "/omega/c\omega   # %1.15f"""%omega+'"  temp1'+label+' > temp2'+label
            cmd += '; cp temp2'+label+' '+'temp1'+label
            #cmd +=  """; sed "/beta/c\beta     # %1.15f"""%betas[j]+'"  temp1'+label+' > temp2'+label
            #cmd += '; cp temp2'+label+' '+'temp1'+label
            cmd +=     """; sed "/mu/c\mu      # %1.15f"""%mu+'" temp1'+label+' > '+input_file_name
            cmd += '; rm temp2'+label+'; rm temp1'+label

            bash_command(cmd)
                
            time.sleep(0.1)
            

            
#save mumap....
