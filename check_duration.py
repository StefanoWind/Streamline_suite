# -*- coding: utf-8 -*-
'''
Check scan durations from hpl files
'''

import os
cd=os.path.dirname(__file__)

import numpy as np
import linecache
import glob 
import pandas as pd
import datetime
from matplotlib import pyplot as plt
import matplotlib.dates as mdates

#%% Inputs
lidar_id=214
root='C:/Lidar/Data/Proc'

date='20230726'

time_res=10#[min]

#%% Initialization
source=os.path.join(date[:4],date[:6],date)

files=glob.glob(os.path.join(root,source,'*.hpl'))

t1=[]
t2=[]
duration=[]
filename=[]
scan_name=[]
date_dt=datetime.datetime.strptime(date, '%Y%m%d')

#%% Main
for f in files:
    try:
        N_tot = sum(1 for line in open(f))
        Nr=int(linecache.getline(f, 3).split(':')[1])
        h1=np.float64(linecache.getline(f, 18).split(' ')[0])
        h2=np.float64(linecache.getline(f, N_tot-Nr).split(' ')[0])
        t1=np.append(t1,date_dt+datetime.timedelta(hours=h1))
        t2=np.append(t2,date_dt+datetime.timedelta(hours=h2)+datetime.timedelta(seconds=1))
        duration=np.append(duration,(h2-h1)*3600)
        scan_name=np.append(scan_name,linecache.getline(f, 8).split(':')[1].split(' - ')[0].strip())
        filename=np.append(filename,f)
        
        linecache.clearcache()
        print(os.path.basename(f)+' done')
    except Exception as e:
        print(e)
    
#%% Ouput
Output=pd.DataFrame()
Output['Filename']=filename
Output['Scan name']=scan_name
Output['Start time']=t1
Output['End time']=t2
Output['Duration [s]']=duration

Output.to_excel('schedules/'+date+str(lidar_id)+'_Schedule.xlsx')


#%% Plots
#scan sequence
sel=scan_name!='Stare'
predefined_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
N=len(np.unique(scan_name[sel]))
cols = predefined_colors[:N]
                                    
plt.figure(figsize=(18,6))

plt.subplot(2,1,1)
t0=np.arange(date_dt,t2[-1],datetime.timedelta(minutes=time_res))
time1=mdates.date2num(t1[sel])
time2=mdates.date2num(t2[sel])
time0=mdates.date2num(t0)

for tt1,tt2,s in zip(time1,time2,scan_name[sel]):
    i_scan=np.where(s==np.unique(scan_name[sel]))[0][0]
    plt.barh(0,tt2-tt1,left=tt1,label=s,color=cols[i_scan])
    
for tt0 in time0:
    plt.plot([tt0,tt0],[-0.5,0.5],'--k')

plt.gca().xaxis_date()
handles, labels = plt.gca().get_legend_handles_labels()
unique_labels = list(set(labels))
unique_handles = [handles[labels.index(label)] for label in unique_labels]
plt.legend(unique_handles, unique_labels)
plt.xlabel('UTC time')
plt.grid()
plt.yticks([0,1],labels=[])
plt.ylim([-0.5,0.5])
plt.xlim(mdates.date2num(t1[0])-5/24/60,mdates.date2num(t2[-1])+5/24/60)

plt.subplot(2,1,2)
sel=scan_name=='Stare'
t0=np.arange(date_dt,t2[-1],datetime.timedelta(minutes=time_res))
time1=mdates.date2num(t1[sel])
time2=mdates.date2num(t2[sel])
time0=mdates.date2num(t0)
ctr=0
for tt1,tt2,s in zip(time1,time2,scan_name[sel]):
    plt.barh(0,tt2-tt1,left=tt1,color='y')
    ctr+=1

for tt0 in time0:
    plt.plot([tt0,tt0],[-0.5,0.5],'--k')

plt.gca().xaxis_date()
plt.xlabel('UTC time')
plt.grid()
plt.yticks([0,1],labels=[])
plt.ylim([-0.5,0.5])
plt.xlim(mdates.date2num(t1[0])-5/24/60,mdates.date2num(t2[-1])+5/24/60)
plt.title('Stare files')
plt.tight_layout()

plt.savefig('schedules/'+date+'_'+str(lidar_id)+'_Schedule.png')


    