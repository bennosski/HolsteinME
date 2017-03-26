import os
import sys

#set up input file for g = sys.argv[1]
g = float(sys.argv[1])
omega = 0.1
N = 64
mu = -N*g**2/omega**2

import subprocess
def bash_command(cmd):
    subprocess.Popen(['/bin/bash', '-c', cmd])

print "setting up input file"

# sed "/chemical/c\mu = $2  #chemical potential" sample_input_forward.txt > "# input_mu_search_$1.txt'

print g

cmd = " cp input temp1"
cmd += """; sed "/g_dqmc/c\%1.15f"""%g+""" # g_dqmc" temp1 > input%1.5f"""%g
print cmd
bash_command(cmd)

print "done with input file"

print "\nsetting up submit file"

mystr = "srun --ntasks=1 --cpus-per-task=1 python ME_holstein.py "

cmd = "cp submit temp3"
cmd += """; sed "/srun/c\\"""+mystr+"input%1.5f"%g+"""" temp3 > submit%1.5f"""%g
print cmd
bash_command(cmd)






