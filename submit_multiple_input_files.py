from numpy import *
import subprocess, time, os, sys


def bash_command(cmd):
    subprocess.Popen(['/bin/bash', '-c', cmd])


submit_ct = 0

mu_map = load('mu_map.npy')
[d1,d2,d3] = shape(mu_map)


def get_count(i,j,k):
    return i*d2*d3 + j*d3 + k

def get_ijk(c):
    i = c/(d2*d3)
    
    r = c%(d2*d3)
    j = r/d3

    k = (c%(d2*d3))%d3

    return i,j,k


start = get_count(0,0,0)
#start = get_count(9,9,9)+1
#end = get_count(5,10,0)
#start = get_count(3,5,6)+1
#end = get_count(0,0,1)
end = get_count(d1-1,d2-1,d3-1)+1

for count in range(start, end):

            [i,j,k] = get_ijk(count)

            
            label = '_%d'%i+'_%d'%j+'_%d'%k
            
            #input_file_name    = 'inputfiles_2_25_17_w4p0_l0p2/input'+label
            #output_folder_name = 'outputfiles_2_25_17_w4p0_l0p2/output'+label+'/'

            p1 = sys.argv[1]
            p2 = sys.argv[2]
            
            #input_file_name    = 'inputfiles_l0p3/input'+label
            #output_folder_name = 'outputfiles_l0p3/output'+label+'/'
            input_file_name = p1+'input'+label
            output_folder_name = p2+'output'+label+'/'
            
            
            #if os.path.exists(output_folder_name+'xsc_Nw40'):
            #    continue

            print 'submitting ',label
            #bash_command('mkdir '+output_folder_name+'; ./multiple_submit '+input_file_name+' '+output_folder_name)
            bash_command('./slac_submit '+output_folder_name)
            #print './slac_submit output_folder_name'

            submit_ct += 1

            #if submit_ct>=882:
            #    1./0
            
            time.sleep(0.1)



print submit_ct
