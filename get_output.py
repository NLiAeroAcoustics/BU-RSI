import shutil
import subprocess
import os
import time
import numpy as np
import pandas as pd
from os.path import exists

skipMissingFile = 1

output_dir = 'examples/rsi_output/ADP_2022/' #output files directory
rsi_sub_dir = 'examples/rsi_run/ADP_2022/' #directory for rsi sub folder

output_name = [] 

for file in os.listdir(rsi_sub_dir):
    d = os.path.join(rsi_sub_dir, file)
    if os.path.isdir(d):
        #print(d)
        output_name.append(d)

#Check files
checklist = []
for i in output_name:
    file_exists = exists(os.path.join(i,'fort.33'))
    checklist.append(file_exists)
    print(i, ' ',file_exists)

if any(checklist) != True :
    print('Warning: Missing files')
    if skipMissingFile == 1:
        quit()

for i in output_name:
    file_exists = exists(os.path.join(i,'fort.33'))
    if file_exists:
        print('loading',i)
        pwl= pd.read_csv(os.path.join(i,'fort.33'), engine ='python',delimiter='        ',names= ['f','pwl_up','pwl_down'])
        pwl = pwl.sort_values(by=['f'])
        outName = i.replace(rsi_sub_dir, '')
        outName = outName.replace('/', '')
        np.savetxt(output_dir+'/'+outName+'.33', pwl.values,fmt='%1.14f')
        print('copied output files for ',i)
    else:
        print('Warning: '+outName+' does not exist!')
    



