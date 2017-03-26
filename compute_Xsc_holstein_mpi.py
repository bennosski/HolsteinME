
from numpy import *
import sys
import os
import time
from init_functions import *

from mpi4py import MPI
comm = MPI.COMM_WORLD
nprocs = comm.size
myrank = comm.rank

t0 = time.time()

folder = sys.argv[1]
Sigma = load(folder+'Sigma.npy')
#the momentum independent Self-energy if this Sigma was calculated using the forward scattering code:
Sigma = Sigma[0,0,:,0,0]

files = os.listdir(folder)
for myfile in files:
    if "input" in myfile:
        inputfile = myfile
        break
        
def parse_line(f):
    line = f.readline()
    index = line.index('#')+1
    if '.' in line:
        return float(line[index:])
    else:
        return int(float(line[index:]))
    
with open(folder+inputfile,'r') as f:
    g_dqmc = parse_line(f)
    Nk     = parse_line(f)
    Nw     = parse_line(f)
    beta   = parse_line(f)
    omega  = parse_line(f)
    superconductivity = parse_line(f)
    mu     = parse_line(f)
    q0     = parse_line(f)
f.close()

g     = g_dqmc * Nk * 1./ sqrt(2. * omega)

if myrank==0:
    print ' g_dqmc ',g_dqmc
    print ' Nk     ',Nk
    print ' Nw     ',Nw
    print ' beta   ',beta
    print ' omega  ',omega
    print ' superconductivty ',superconductivity
    print ' q0     ',q0
    print ' mu     ',mu

q0    = 2*pi*q0

NwOriginal = Nw
Nws = [20,30,40,50]

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
         print 'done with x0 \n'

     if abs(g) > 1e-10:  # g is nonzero         
         g2D = zeros([Nw,Nw], dtype=complex)
         for n1 in range(Nw):
             for n2 in range(Nw):
                 g2D[n1,n2] = -2.*g**2*omega / ((wn[n1]-wn[n2])**2+omega**2)

         t = g2D.copy()
         change = 1.0

         iter_selfconsistency = 20
         for myiter in range(iter_selfconsistency):
             if change < 1e-4:
                 break

             if myrank==0:
                 print 'iter ',myiter

             tnew = zeros([Nw,Nw], dtype=complex)
             tnew_proc = zeros([Nw,Nw], dtype=complex)

             for ik1 in range(Nk):
                 for ik2 in range(Nk):
                     if (ik1+ik2*Nk)%nprocs == myrank:
                         for n in range(Nw):
                             for np in range(Nw):
                                 for npp in range(Nw):
                                     tnew_proc[np,n] -= 1./(beta*Nk**2)*1./(Z[npp]**2*wn[npp]**2 + band[ik1,ik2]**2)*g2D[npp,np] * t[npp,n]

             comm.Allreduce(tnew_proc, tnew, op=MPI.SUM)

             tnew += g2D 
             change = sum(abs(tnew-t))/Nw**2
             if myrank==0:
                 print ' '
                 print change
             t = tnew.copy()

             if myrank==0:
                 print '\ntime elapsed ',time.time()-t0

         #save(folder,'t')

         if myrank==0:
             print 'computing xsc'

         x = asarray(0., dtype=complex)
         x_proc = asarray(0., dtype=complex)
         for ik1 in range(Nk):
             for ik2 in range(Nk):
                 if (ik1 + ik2*Nk)%nprocs==0:
                     for ip1 in range(Nk):
                         for ip2 in range(Nk):
                             for n in range(Nw):
                                 for np in range(Nw):
                                     x_proc -= 1./(beta*Nk**2)**2 * 1./(Z[np]**2*wn[np]**2 + band[ip1,ip2]**2) * t[np, n] * 1./(Z[n]**2*wn[n]**2 + band[ik1,ik2]**2)  

         comm.Allreduce(x_proc, x, op=MPI.SUM)

         x += x0
     else: # only the bare part
         x = x0
         
     if myrank==0:
         print '---------------------'
         print 'Xsc = ', x
         print '---------------------'

         print 'total time elapsed ',time.time()-t0

         savetxt(folder+'xsc_Nw%d'%Nw, [real(x)])
    
