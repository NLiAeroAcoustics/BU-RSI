############RSI python wrapper##########
####Import matlab, pgi before running this code####
import shutil
import subprocess
import os
import time
import numpy as np

###########Inputs#################

parafreq = 3 #run multiple frequencies at a time ? 0=no,otherwise number of processes to run at a time
max_runtime = '140:00:00' #maximum run time, only used if parafreq = 0
#Adjust frequency spectrum when needed. number of frequencies between
#1-800Hz, 800-2000Hz, 2000-3000Hz, 3000-4000Hz, 4000-8000Hz,
#8000-12000Hz, 12000-16000Hz, 16000-21000Hz, 21000-25000Hz
nfreq = [8,5,4,3,2,2,2,4,3]

nf = 36 #Number of frequencies to be computed only used if parafreq = 0
input_dir = 'examples/rsi_input/ADP_2022/' #input files directory
findfreq_dir = 'findfreqs_src/' #frequency selector directory
rsi_main_dir = 'rsi_src' #directory for rsi source code
rsi_sub_dir = 'examples/rsi_run/ADP_2022/' #directory for rsi run folder
# Inputs below are only available when parafreq > 1
maxf = 20000 #Maximum frequency
ng = 300 #Number of k3(spanwise) integration points
nc = 80  #Number of chordwise integration points
nrcalc = 45 #Number of radial strips
input_name = [] #(Optional) Specify case names. Leave blank to run all cases 

#################################
## If input_name is not specified, run all cases in the directory
if not input_name:
    for file in os.listdir(input_dir):
        d = os.path.join(input_dir, file)
        if 'aero95const.input' in d:
            #print(d)
            d = d.replace(input_dir, '')
            d = d.replace('/', '')
            d = d.replace('_aero95const.input', '')
            input_name.append(d)

    input_name = np.unique(input_name)
print('number of cases =',len(input_name))
parent_dir = os.getcwd()
#os.system('module load matlab')
for i in input_name:
    isExist = os.path.exists(rsi_sub_dir+'/'+i)
    if not isExist:
        shutil.copytree(rsi_main_dir,rsi_sub_dir +'/'+str(i)) 
        shutil.copy(os.path.join(input_dir,i+'_aero95const.input'),os.path.join(rsi_sub_dir,i,'aero95const.input'))
        shutil.copy(os.path.join(input_dir,i+'_turb95mat.input'),os.path.join(rsi_sub_dir,i,'turb95mat.input'))
        shutil.copy(os.path.join(input_dir,i+'_aero95mat.input'),os.path.join(rsi_sub_dir,i,'aero95mat.input'))
        print('copied input files for ',i)
        with open(os.path.join(rsi_sub_dir,i,'aero95const.input'), "r") as aero:
            val = aero.readlines()
            line1 = val[1].split(",")
            nv = line1[0].split("=")
            nv = nv[1]
            nb = line1[1].split("=")
            nb = nb[1]
            nrdat = line1[2].split("=")
            nrdat = nrdat[1]
            
            t = val[3].split(",")
            t = t[0].split("=")
            t = t[1]
            m_ann = val[3].split(",")
            m_ann = m_ann[2].split("=")
            m_ann = m_ann[1]
        print(t,m_ann)
        with open(os.path.join(findfreq_dir,'freqs.m'), "r") as findfreqs:
            line = findfreqs.readlines()
            line[3] = 'M_ann =' + str(m_ann)+'; \n'
            line[5] = 'T =' + str(t)+'; \n'
            line[12] = 'nfreq0 =' + str(nfreq[0])+'; \n'
            line[13] = 'nfreq1 =' + str(nfreq[1])+'; \n'
            line[14] = 'nfreq2 =' + str(nfreq[2])+'; \n'
            line[15] = 'nfreq3 =' + str(nfreq[3])+'; \n'
            line[16] = 'nfreq4 =' + str(nfreq[4])+'; \n'
            line[17] = 'nfreq5 =' + str(nfreq[5])+'; \n'
            line[18] = 'nfreq6 =' + str(nfreq[6])+'; \n'
            line[19] = 'nfreq7 =' + str(nfreq[7])+'; \n'
            line[20] = 'nfreq8 =' + str(nfreq[8])+'; \n'
        with open(os.path.join(findfreq_dir,'freqs2.m'), "w") as findfreqs:
            findfreqs.writelines(line)
        print('running freqs2.m for ',i)
        os.system('''matlab -nodisplay -nosplash -nodesktop -r "run('''+"'"+ findfreq_dir +'''/freqs2.m');exit;"''')
        os.system('module load pgi')
        shutil.copy(os.path.join(findfreq_dir,'fs.dat'),os.path.join(rsi_sub_dir,i,'fs.dat'))
        print('copied fs.dat to rsi_sub for ',i)
        os.chdir(os.path.join(parent_dir,rsi_sub_dir,i))
        with open("aero95const.input", "r") as f:
            lines = f.readlines()
        lines[1] = "NV = "+str(nv)+", NB = "+str(nb)+", NRdat="+str(nrdat)+", NRcalc = "+str(nrcalc)+", NC = "+str(nc)+","+"\n"
        with open("aero95const.input", "w") as f:
            f.writelines(lines)        
        if parafreq >= 1:
            with open('batchrsi.py', "r") as batchrsi:
                line = batchrsi.readlines()
                line[5] = 'maxf =' +str(maxf) +' \n'
                line[6] = 'core =' + str(parafreq)+' \n'
                line[7] = 'ng =' + str(ng)+' \n'
            with open('batchrsi.py', "w") as batchrsi:
                batchrsi.writelines(line)
            os.system('python3 batchrsi.py')
        if parafreq == 0:
            with open("aero95const.input", "r") as f:
                lines = f.readlines()
                lines[2] = 'NF = ' + str(nf) + ', Ng = 200,' +"\n"
            with open("aero95const.input", "w") as f:
                f.writelines(lines)
            os.system('qsub -b y -l h_rt='+max_runtime+' rsi95')
        os.chdir(parent_dir)
    else:
        print(str(rsi_sub_dir+'/'+i)+' already exist!')
