# -*- coding: utf-8 -*-
#03/28/2023: created
#03/31/2032: adapted to new Halo simulator, finalized
#04/04/2023: added volumetric option
#07/05/2023 (v 2): fixed bug when azi and ele have different lengths

import os
cd=os.path.dirname(__file__)
from datetime import datetime
import utils as utl
import Halo_scan_sim as HS
import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt

plt.close('all')
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm' 
matplotlib.rcParams['font.size'] = 14

#%% Inputs
source='scans/250124_awaken_profiling/mod_inc_six-beam.xlsx'
source_time='data/Halo_time_info.xlsx'
volumetric=False
scan_mode='SSM'
identifier='no-overlapping'
model='XR+'
ppr=10000
repeat=1

#%% Initialization
scan_info=pd.read_excel(source)
azi=scan_info['Azimuth'].values
ele=scan_info['Elevation'].values

if volumetric:
    [azi,ele]=np.meshgrid(azi,ele)
    azi=azi.ravel()
    ele=ele.ravel()
    
azi=np.append(azi,azi[0])
ele=np.append(ele,ele[0])

sel=~np.isnan(azi+ele)
azi=azi[sel]
ele=ele[sel]
    
L=''
r=np.arange(0,1,0.1)

name='scans/'+datetime.strftime(datetime.now(),'%y%m%d_%H%M')+'_'+os.path.basename(source[:-5])+'_'+model+'_'+scan_mode+'_'+identifier+'_'+str(ppr)+'x'+str(repeat)

#%% Main

#read kinematic parameters
time_info=pd.read_excel(source_time,sheet_name=model)
time_info=time_info.set_index('Scan settings')

tp=time_info['Processing time'][scan_mode+ ' - '+identifier]
ta=time_info['Acquisition time'][scan_mode+ ' - '+identifier]
tau_s=tp+ta*ppr
time=[0]

if scan_mode=='CSM':
    
    ppd1=500000/360#points per degree
    ppd2=250000/360#points per degree
    
    dazi=np.diff(azi)
    dele=np.diff(ele)
    stop=np.concatenate(([0],np.where(np.abs((np.diff(dazi)+np.diff(dele)))>10**-10)[0]+1,[-1]))
    
    th_range=azi[stop]
    P1=-th_range*ppd1
    azi2=[th_range[0]]
    s1=np.abs(dazi[stop[:-1]]/tau_s)
    s1[s1>50000/ppd1]=50000/ppd1
    S1=np.append(5000,s1*ppd1/10)
    S1[S1>5000]=5000
    A1=50+np.zeros(len(P1))
    a1=A1*1000/ppd1
    
    b_range=ele[stop]
    P2=-b_range*ppd2
    ele2=[b_range[0]]
    s2=np.abs(dele[stop[:-1]])/tau_s
    s2[s2>50000/ppd2]=50000/ppd2
    S2=np.append(5000,s2*ppd2/10)
    S2[S2>5000]=5000
    A2=50+np.zeros(len(P2))
    a2=A2*1000/ppd2

x,y,z=utl.lidar_xyz(r, ele, azi)

#write scan file

if scan_mode=='CSM':
    for p1,p2,s1,s2,a1,a2 in zip(P1[:-1],P2[:-1],S1[:-1],S2[:-1],A1[:-1],A2[:-1]):
        
        l='A.1=%.0f'%a1+',S.1=%.0f'%s1+',P.1=%.0f'%p1+'*A.2=%.0f'%a2+',S.2=%.0f'%s2+',P.2=%.0f'%p2+'\nW=0\n'
        L+=l
elif scan_mode=='SSM':
    if utl.len2(azi)>1:
        for a,e in zip(azi[:-1],ele[:-1]):
            L=L+('%07.3f' % a+ '%07.3f' % e +'\n')
    else:
        L='%07.3f' % azi+ '%07.3f' % ele +'\n'


with open(name+'.txt','w') as fid:
    fid.write(L*repeat)
    fid.close()
Output=pd.DataFrame()


#simulate scanning head
time_sim,azi_sim,ele_sim=HS.Halo_scan_sim(name+'.txt',ppr,identifier,model,source_time,dt=0.01,ang_tol=0.25,azi0=0,ele0=0, ppd1=500000/360,ppd2=250000/360)

Output['Time']=time
Output['Azimuth']=azi2
Output['Elevation']=ele2
Output.set_index('Time')
Output.to_excel(name+'.xlsx')

#%% Plots

plt.figure(figsize=(18,8))
plt.subplot(2,1,1)
plt.plot(time,azi2,'.-k')
plt.xlabel('Time [s]')
plt.ylabel(r'$\theta$ [$^\circ]$')
plt.grid()
plt.subplot(2,1,2)
plt.plot(time,ele2,'.-k')
plt.xlabel('Time [s]')
plt.ylabel(r'$\beta$ [$^\circ]$')
plt.grid()

fig=plt.figure()
ax = plt.subplot(1,1,1,projection='3d')
ax.scatter(x,y,z,s=5,c='k', alpha=0.25)
ax.set_xlabel(r'$x$ [m]')
ax.set_ylabel(r'$y$ [m]')
ax.set_zlabel(r'$z$ [m]')
utl.axis_equal()
plt.grid()