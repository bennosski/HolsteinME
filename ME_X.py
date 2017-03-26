# -*- coding: utf-8 -*-
"""
Created on Thu Dec 08 01:44:04 2016

@author: Ben
"""

import subprocess
def bash_command(cmd):
    subprocess.Popen(['/bin/bash', '-c', cmd])

#bash_command('source mlpython2.7.5.sh')
    
import numpy as np
from numpy import *
from init_functions import *
from Functions import *
import time
import subprocess
import sys

myrank = 0
nprocs = 1

###Nk must be odd or else the momentum points do not 
###form a group under addition!
###Choose Nw to be even (number of Fermionic frequencies)


# argv[1] is the input file
# argv[2] is the folder where everything is saved


# g_ME = g_dqmc * sqrt(N_total_DQMC) / sqrt(2 omega)

#bash_command('source setup.sh')

tstart = time.time()

def parse_line(f):
    line = f.readline()
    index = line.index('#')+1
    if '.' in line:
        return float(line[index:])
    else:
        return int(float(line[index:]))
    
with open(sys.argv[1],'r') as f:
    g_dqmc = parse_line(f)
    Nk     = parse_line(f)
    Nw     = parse_line(f)
    beta   = parse_line(f)
    omega  = parse_line(f)
    superconductivity = parse_line(f)
    mu = parse_line(f)
    q0     = parse_line(f)
f.close()

#savedir = sys.argv[2]
import os
'''
from datetime import date
today = date.today()
yr = today.timetuple()[0]
mn = today.timetuple()[1]
dy = today.timetuple()[2]
savedir = 'data_%d'%mn+'_%d'%dy+'_%d'%yr+'/fig3/q%1.1f'%q0+'_omega%1.1f'%omega+'_g%1.3f'%g_dqmc+'_mu%1.3f'%abs(mu)+'_Nw%d'%Nw+'_Nk%d'%Nk + '_beta%1.1f'%beta+'/'
'''
savedir = sys.argv[2]
print 'savedir = ',savedir
if not os.path.exists(savedir) and myrank==0:
    os.makedirs(savedir)

if myrank==0:
    bash_command('cp '+sys.argv[1]+' '+savedir)

    print ' g_dqmc ',g_dqmc
    print ' Nk     ',Nk
    print ' Nw     ',Nw
    print ' beta   ',beta
    print ' omega  ',omega
    print ' superconductivty ',superconductivity
    print ' mu     ',mu
    print ' q0     ',q0

#comm.Barrier()

print "hi ",myrank

g     = g_dqmc * Nk * 1./ sqrt(2. * omega)

q0    = 2*pi*q0

iter_selfconsistency = 300

kxs, kys  = init_momenta(Nk)
gofq      = init_gofq(kxs, kys, Nk, g, q0)
iw_bose   = init_boson_freq(Nw, beta)
iw_fermi  = init_fermion_freq(Nw, beta)
band      = init_band(kxs, kys, Nk, mu)
D         = init_D(Nw, beta, omega, iw_bose)

#now do the same calculation but with the FFT
#G         = init_G(Nk, Nw, beta, omega, band, kxs, kys, iw_fermi, superconductivity)
#G         = load("data/G.npy")

G         = zeros([Nk,Nk,Nw,2,2], dtype=complex)
G_proc    = zeros([Nk,Nk,Nw,2,2], dtype=complex)
Conv      = zeros([Nk,Nk,2,2], dtype=complex)
Sigma     = zeros([Nk,Nk,Nw,2,2], dtype=complex)

fft_gofq2 = fft.fft2(gofq**2)

change = ones([2,2], dtype=complex)

if superconductivity:
    Sigma[:,:,:,0,1] = 0.01*1j
    Sigma[:,:,:,1,0] = 0.01*1j


if os.path.exists(savedir+'Sigma.npy'):
    Sigma = load(savedir+'Sigma.npy')
    iter_selfconsistency = 0
    print 'loading Sigma from previous run. Not performing selfconsistency'
    
#selfconsistency loop
for myiter in range(iter_selfconsistency):
    if abs(change[0,0]) < 1e-10:
        break
    
    #compute new G
    for ik1 in range(Nk):
        for ik2 in range(Nk):
            
            for n in range(Nw):
                iwn = iw_fermi[n]
                
                G_proc[ik1,ik2,n,:,:] = linalg.inv(iwn*tau0 - band[ik1,ik2]*tau3 - Sigma[ik1,ik2,n,:,:])

    G = zeros([Nk,Nk,Nw,2,2], dtype=complex)
    G = G_proc
    #comm.Allreduce(G_proc, G, op=MPI.SUM)

    
    Sigma_old = Sigma.copy()
    Sigma_proc = zeros([Nk,Nk,Nw,2,2], dtype=complex)
    Sigma = zeros([Nk,Nk,Nw,2,2], dtype=complex)
    
    #compute new Sigma    
    change = zeros([2,2], dtype=complex)

    for m in range(Nw): ### check this bound. it was n
        fft_G = fft.fft2(einsum('ij,abjk,kl->abil',tau3,G[:,:,m,:,:],tau3), axes=(0,1))

        for n in range(Nw): 
            n_m = subtract_freqs(n,m,Nw)
                            
            if(n_m>=0 and n_m<Nw-1):

                Conv = 1.0/(Nk**2)/beta * D[n_m] * fft.ifft2( einsum('ij,ijab->ijab', fft_gofq2 , fft_G) , axes=(0,1))
                Conv = roll(Conv, -Nk/2, axis=0)
                Conv = roll(Conv, -Nk/2, axis=1)

                Sigma_proc[:,:,n,:,:] -= Conv
                                    
    Sigma = Sigma_proc
    #comm.Allreduce(Sigma_proc, Sigma, op=MPI.SUM)
    
    change += einsum('ijklm->lm',abs(Sigma-Sigma_old))/Nk**2
    dens = abs(1.0 + 2.0*einsum('ijk->',G[:,:,:,0,0])/Nk**2/beta)
    
    if myrank==0:
        print " "
        print "iteration ",myiter
        print "change ", change
        print "iteration time ",time.time() - tstart
        print "filling : ", dens

        
if myrank==0:
    save(savedir+"GM.npy", G)
    save(savedir+"Sigma.npy", Sigma)    
    savetxt(savedir+"dens",[dens])
    print "---------------------------"
    print "total run time ", time.time() - tstart
    print "Done with Migdal piece"
    print "Starting X calculation"
    print " "

    
# copy input file into the savedir

t0 = time.time()

#folder = sys.argv[1]
#Sigma = load(folder+'Sigma.npy')
#the momentum independent Self-energy if this Sigma was calculated using the forward scattering code:
Sigma = Sigma[0,0,:,0,0]


NwOriginal = Nw
#Nws = [20,30,40,50]
Nws = [20, 30, 40]

for Nw in Nws:

     kxs, kys = init_momenta(Nk)
     band     = init_band(kxs, kys, Nk, mu)
     iwn      = init_fermion_freq(Nw, beta)
     wn = imag(iwn)

     mid = NwOriginal/2
     Z = 1.0 - Sigma[mid-Nw/2:mid+Nw/2]/iwn

     x0 = 0.
     for ik1 in range(Nk):
         for ik2 in range(Nk):
             for n in range(Nw):
                 x0 += 1./(beta*Nk**2) * 1./(Z[n]**2*wn[n]**2 + band[ik1,ik2]**2)

     if myrank==0:
         print 'done with x0 ', x0

     if abs(g) > 1e-10:  # g is nonzero         
         g2D = zeros([Nw,Nw], dtype=complex)
         for n1 in range(Nw):
             for n2 in range(Nw):
                 g2D[n1,n2] = -2.*g**2*omega / ((wn[n1]-wn[n2])**2+omega**2)

         t = g2D.copy()
         change = 1.0

         iter_selfconsistency = 20
         for myiter in range(iter_selfconsistency):
             if change < 1e-6:
                 break

             if myrank==0:
                 print 'iter ',myiter

             tnew = zeros([Nw,Nw], dtype=complex)
             tnew_proc = zeros([Nw,Nw], dtype=complex)

             for ik1 in range(Nk):
                 for ik2 in range(Nk):
                     for n in range(Nw):
                         for np in range(Nw):
                             for npp in range(Nw):
                                 tnew_proc[np,n] -= 1./(beta*Nk**2)*1./(Z[npp]**2*wn[npp]**2 + band[ik1,ik2]**2)*g2D[npp,np] * t[npp,n]

             #comm.Allreduce(tnew_proc, tnew, op=MPI.SUM)
             tnew = tnew_proc
             tnew += g2D 
             change = sum(abs(tnew-t))/Nw**2 / sum(abs(tnew))
             if myrank==0:
                 print 'change: ',change
                 print 't sum: ',sum(abs(tnew))
             t = tnew.copy()

             if myrank==0:
                 print 'iter time elapsed ',time.time()-t0
                 print ' '
                 
         #save(savedir+'t', t)

         if myrank==0:
             print 'computing xsc'

         x = asarray(0.)
         x_proc = asarray(0.)
         for ik1 in range(Nk):
             for ik2 in range(Nk):
                 for ip1 in range(Nk):
                     for ip2 in range(Nk):
                         for n in range(Nw):
                             for np in range(Nw):
                                 x_proc -= real( 1./(beta*Nk**2)**2 * 1./(Z[np]**2*wn[np]**2 + band[ip1,ip2]**2) * t[np, n] * 1./(Z[n]**2*wn[n]**2 + band[ik1,ik2]**2) )

         #comm.Allreduce(x_proc, x, op=MPI.SUM)
         x = x_proc

         print 'xs'
         print x, x0
         
         x += x0
     else: # only the bare part
         x = x0
         
     if myrank==0:
         print '---------------------'
         print 'Xsc = ', x
         print '---------------------'

         print 'total time elapsed ',time.time()-t0

         savetxt(savedir+'xsc_Nw%d'%Nw, [real(x)])

