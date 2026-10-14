import shutil
import subprocess
import os
import time
import numpy as np
maxf = 20000
core =16
ng = 200
n = list(range(1,core+1))
command = ["rsi95"]
freq = []
freq_residual = []
with open("fs.dat", 'r') as out_file:
    freqlines = out_file.readlines()
for j in range(len(freqlines)):
    slice_  = freqlines[j]
    splitslice = slice_.split()
    freq.append(float(splitslice[0]))
for i in range(len(freq)):
    if float(freq[i]) >= maxf:
        del freq[i+1:-1]
        break
if freq[-1] and freq[-2] >= maxf:
    del freq[-1]
print('Computing frequencies, ', freq)
residual = int(len(freq) % core)    
for i in range(0,core):    
    if len(freq) % core != 0:
        freq_residual.append(freq[-1])
        del freq[-1]
    i + 1

freq = np.array(freq).reshape(int(len(freq)/core),core).transpose()
shutil.copyfile("fs.dat", "fs0.dat")
nfreq = len(freq[1,:])
with open("aero95const.input", "r") as f:
        lines = f.readlines()
        lines[2] = "NF = 0, Ng = "+str(ng)+ ", \n"
with open("aero95const.input", "w") as f:
    f.writelines(lines)
            

for i in range(0,core):
    shutil.copyfile("fs.dat", "fs"+str(i)+".dat")
    with open("fs"+str(i)+".dat", "r") as f:
        lines = f.readlines()
    with open("fs"+str(i)+".dat", "w") as f:
        for k in freq[i,:]:
            f.write(str(k)+"\n")
    f.close()
    
if residual != 0:
    for i in range(residual):
        with open("fs"+str(i)+".dat", "r") as f:
            lines = f.readlines()
            lines.append(str(freq_residual[i]))
            #print(lines)
        with open("fs"+str(i)+".dat", "w") as f:
            f.writelines(lines)
            
for i in range(0,core):    
    with open('rsi.f95', 'r') as rsi95:
        rsi_lines = rsi95.readlines()
        rsi_lines[56] ='''open (unit = 22, file = 'fs''' +str(i)+'''.dat', form = 'formatted',  &\n'''
    with open('rsi.f95', 'w') as rsi95:
        rsi95.writelines( rsi_lines )
    os.system("make -f Makersifull.pgi")
    time.sleep(5)
    shutil.move("rsi95", "rsi"+str(i))
    os.system("qsub -b y -l cpu_arch='cascadelake|skylake|icelake|sapphirerapids' -l h_rt=12:00:00 rsi"+str(i))
    #os.system("qsub -b y -l cpu_arch='cascadelake|skylake|icelake|sapphirerapids' rsi"+str(i))
    print('process',i+1,'submitted')
    
