
import subprocess
def bash_command(cmd):
    subprocess.Popen(['/bin/bash', '-c', cmd])

import os
import sys

dirpath = sys.argv[1]
folders = os.listdir(dirpath)

for folder in folders:
    print folder
    newf = folder + '_for_download/'
    cmd  = 'cp -r '+dirpath+folder+' '+dirpath+newf
    cmd += '; rm '+dirpath+newf+'GM.npy'
    cmd += '; rm '+dirpath+newf+'Sigma.npy'
    bash_command(cmd)
    
    
