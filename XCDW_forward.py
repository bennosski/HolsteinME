# -*- coding: utf-8 -*-
"""
Created on Thu Dec 08 01:44:04 2016

@author: Ben
"""

import subprocess
def bash_command(cmd):
    subprocess.Popen(['/bin/bash', '-c', cmd])

#bash_command('source mlpython2.7.5.sh')
    
from numpy import *
from init_functions import *
from Functions import *
import time
import subprocess
import sys
import os

myrank = 0
nprocs = 1


savedir = sys.argv[1]
print 'savedir = ',savedir


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

myfiles = os.listdir(savedir)
for myfile in myfiles:
    if 'input' in myfile:
        with open(savedir+myfile,'r') as f:
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
'''
from datetime import date
today = date.today()
yr = today.timetuple()[0]
mn = today.timetuple()[1]
dy = today.timetuple()[2]
savedir = 'data_%d'%mn+'_%d'%dy+'_%d'%yr+'/fig3/q%1.1f'%q0+'_omega%1.1f'%omega+'_g%1.3f'%g_dqmc+'_mu%1.3f'%abs(mu)+'_Nw%d'%Nw+'_Nk%d'%Nk + '_beta%1.1f'%beta+'/'
'''

if myrank==0:
    print ' g_dqmc ',g_dqmc
    print ' Nk     ',Nk
    print ' Nw     ',Nw
    print ' beta   ',beta
    print ' omega  ',omega
    print ' superconductivty ',superconductivity
    print ' mu     ',mu
    print ' q0     ',q0

#comm.Barrier()


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


gofq2 = gofq**2
#fft_gofq2 = fft.fft2(gofq**2)


'''
if myrank==0:
    save(savedir+"GM.npy", G)
    save(savedir+"Sigma.npy", Sigma)    
    savetxt(savedir+"dens",[dens])
    print "---------------------------"
    print "total run time ", time.time() - tstart
    print "Done with Migdal piece"
    print "Starting X calculation"
    print " "
'''

Sigma = load(savedir+"Sigma.npy")
     
# copy input file into the savedir

t0 = time.time()

#folder = sys.argv[1]
#Sigma = load(folder+'Sigma.npy')
#the momentum independent Self-energy if this Sigma was calculated using the forward scattering code:
Sigma = Sigma[:,:,:,0,0]

NwOriginal = Nw
#Nws = [20,30,40,50]
Nws = [20, 30, 40]

for Nw in Nws:

    Xcdw = zeros([5,5])

    for iQ1 in range(4,-1,-1):
        for iQ2 in range(iQ1,-1,-1):
            kxs, kys = init_momenta(Nk)
            band     = init_band(kxs, kys, Nk, mu)
            iwn      = init_fermion_freq(Nw, beta)
            wn = imag(iwn)

            mid = NwOriginal/2
            Z = 1.0 - Sigma[:,:,mid-Nw/2:mid+Nw/2]/iwn

            bubble = zeros([Nw,Nk,Nk],dtype=complex)
            for n in range(Nw):
                for ik1 in range(Nk):
                    for ik2 in range(Nk):
                        kpq1 = (ik1+iQ1)%Nk
                        kpq2 = (ik2+iQ2)%Nk                        
                        #bubble[n,ik1,ik2] = 1./(Z[ik1,ik2,n]*wn[n] - band[ik1,ik2]) * 1./(Z[kpq1,kpq2,n]*wn[n] + band[kpq1,kpq2])
                        bubble[n,ik1,ik2] = -1./(Z[ik1,ik2,n]*wn[n] - band[ik1,ik2]) * 1./(Z[kpq1,kpq2,n]*wn[n] - band[kpq1,kpq2])


            print 'made bubble'
                        
            x0 = 0.
            for ik1 in range(Nk):
                for ik2 in range(Nk):
                    for n in range(Nw):
                        #x0 += 1./(beta*Nk**2) * 1./(Z[n]**2*wn[n]**2 + band[ik1,ik2]**2)
                        x0 += 1./(beta*Nk**2) * bubble[n,ik1,ik2]

            if myrank==0:
                print 'done with x0 ', x0

            if abs(g) > 1e-10:  # g is nonzero         
                g2D = zeros([Nk,Nk,Nw,Nw], dtype=complex)
                for iq1 in range(Nk):
                    for iq2 in range(Nk):
                        for n1 in range(Nw):
                            for n2 in range(Nw):
                                g2D[iq1,iq2,n1,n2] = -2.*gofq2[iq1,iq2]*omega / ((wn[n1]-wn[n2])**2+omega**2)

                t = g2D.copy()
                change = 1.0

                iter_selfconsistency = 20
                for myiter in range(iter_selfconsistency):
                    if change < 1e-6:
                        break

                    if myrank==0:
                        print 'iter ',myiter

                    tnew = zeros([Nk,Nk,Nw,Nw], dtype=complex)
                    tnew_proc = zeros([Nk,Nk,Nw,Nw], dtype=complex)

                    for iq1 in range(Nk):
                        for iq2 in range(Nk):
                            for ik1 in range(Nk):
                                for ik2 in range(Nk):
                                    for n in range(Nw):
                                        for np in range(Nw):
                                            for npp in range(Nw):
                                                #tnew_proc[np,n] -= 1./(beta*Nk**2)*1./(Z[npp]**2*wn[npp]**2 + band[ik1,ik2]**2)*g2D[npp,np] * t[npp,n]
                                                tnew_proc[iq1,iq2,np,n] -= 1./(beta*Nk**2) * bubble[npp,ik1,ik2] * g2D[iq1,iq2,npp,np] * t[iq1,iq2,npp,n]

                    #comm.Allreduce(tnew_proc, tnew, op=MPI.SUM)
                    tnew = tnew_proc
                    tnew += g2D 
                    change = sum(abs(tnew-t))/Nw**2 / sum(abs(tnew))
                    if myrank==0:
                        print change
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
                                iqx, iqy = subtract_momenta(ip1,ip2,ik1,ik2,Nk)
                                for n in range(Nw):
                                    for np in range(Nw):
                                        #x_proc -= real( 1./(beta*Nk**2)**2 * 1./(Z[np]**2*wn[np]**2 + band[ip1,ip2]**2) * t[np, n] * 1./(Z[n]**2*wn[n]**2 + band[ik1,ik2]**2) )
                                        x_proc -= real( 1./(beta*Nk**2)**2 * bubble[np,ip1,ip2] * t[iqx,iqy,np, n] * bubble[n,ik1,ik2] )

                #comm.Allreduce(x_proc, x, op=MPI.SUM)
                x = x_proc

                print 'xs'
                print x, x0

                x += x0
            else: # only the bare part
                x = x0

            if myrank==0:
                print '---------------------'
                print 'Xcdw = ', x
                print '---------------------'

                print 'total time elapsed ',time.time()-t0

                Xcdw[iQ1,iQ2] = real(x)
                
                
    savetxt(savedir+'xcdw_Nw%d'%Nw, Xcdw)
        
