'''
Build daily scan schedule
'''

import os
cd=os.path.dirname(__file__)
import datetime as dt
import numpy as np
import pandas as pd

#%% Inputs
source='scans/250509_awaken_multiPPT'
source_scan='seq_awaken_multiPPT.xlsx'
time_corr=-2
avg_loops=1
scan_mode='S'

#%% Initialization
scan_info=pd.read_excel(os.path.join(cd,source,source_scan)) 

filenames=scan_info['Filename'].values
T=scan_info['Total time [s]'].values
dT=scan_info['Sampling time [s]'].values

N=np.floor(T*(1+time_corr/100)/dT).astype(int)

scan_info['Reps']=N
scan_info['Total time (real)']=dT*N

azi_all=[]
ele_all=[]
Nb=[]

now=dt.datetime.now()
dir_save='scans/'+dt.datetime.strftime(now,'%y%m%d_%H%M')+'_'+os.path.basename(source_scan)[:-5]

#%% Main
os.makedirs(dir_save,exist_ok=True)

for f,n in zip(filenames,N):
    with open(os.path.join(source,f),'r') as fid:
        lines=fid.read()
    
    s=lines*n
    
    with open(os.path.join(dir_save,os.path.basename(f[:-4])+'x'+str(n)+'.txt'),'w') as fid:
        fid.write(s[:-1])
        

N_day=np.floor(24*3600/np.sum(T)).astype(int)
s=''
for i_seq in range(N_day):
    for i_scan in range(len(T)):
        t=np.float64(i_seq*np.sum(T)+np.sum(T[:i_scan]))
        
        t_str=dt.datetime.strftime(dt.datetime(2000, 1, 1)+dt.timedelta(seconds=t),'%H%M%S')
        f=os.path.basename(filenames[i_scan][:-4])+'x'+str(N[i_scan])
        s+=t_str+'\t'+f+'\t'+str(avg_loops)+'\t'+scan_mode+'\t7\n'
        
with open(dir_save+'.dss','w') as fid:
    fid.write(s[:-1])
    
    
scan_info.set_index('Scan name').to_excel(source+source_scan)