# -*- coding: utf-8 -*-
'''
Linear model of accumulation time
'''

import os
cd=os.path.dirname(__file__)
import sys
import yaml
import numpy as np
import linecache
import glob 
from matplotlib import pyplot as plt
import matplotlib
plt.close('all')
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm' 
matplotlib.rcParams['font.size'] = 14

#%% Inputs
source_config='config.yaml'
source='CSM_test/*hpl'#data folder
percentile=50#percentile of the sampling time selected as the nominal one
ang_tol=0.25#[deg] angular tolerance

#graphics
max_time=450#[s] maximum time in the plot

#%% Initialization
with open(source_config, 'r') as fid:
    config = yaml.safe_load(fid)
    
#imports
sys.path.append(config['path_utils'])
import utils as utl
files=glob.glob(os.path.join(cd,'data',source))
tau_nom=[]
tau_min=[]
tau_max=[]
ppr=[]

#%% Main
for f in files:      
    t=[]
    ele=[]
    azi=[]
    ctr2=0
    
    #read hpl file
    N_tot = sum(1 for line in open(f))
    Nr=int(linecache.getline(f, 3).split(':')[1])
    ppr=np.append(ppr,int(linecache.getline(f, 6).split(':')[1]))
    
    line_number=18
    while line_number<N_tot:
        l=linecache.getline(f, line_number).split(' ')
        l_strip=[]
        for ll in l:
            if len(ll)>0:
                l_strip.append(ll)
        
        t=np.append(t,np.float64(l_strip[0])*3600) 
        azi=np.append(azi,np.float64(l_strip[1]))
        ele=np.append(ele,np.float64(l_strip[2]))          
        ctr2+=1
        line_number=18+ctr2*(Nr+1)          
    t=t-t[0]
    
    #accumulation time estimation
    tau_nom=np.append(tau_nom,np.nanpercentile(np.diff(t),percentile))#nominal 
    tau_min=np.append(tau_min,np.nanpercentile(np.diff(t),5))#max 
    tau_max=np.append(tau_max,np.nanpercentile(np.diff(t),95))#min 
    linecache.clearcache()
    
    #plot trajectory
    azi=utl.round(azi,ang_tol)%360
    ele=utl.round(ele,ang_tol)%360
    plt.figure(figsize=(18,8))
    plt.subplot(2,1,1)
    plt.plot(t,azi,'.-k',markersize=1,linewidth=1)
    plt.ylabel(r'Azimuth [$^\circ$]')
    plt.xlim([0,max_time])
    plt.grid()
    plt.title(os.path.basename(f)+', PPR='+str(int(ppr[-1])))
    
    plt.subplot(2,1,2)
    plt.plot(t,ele,'.-k',markersize=1,linewidth=1)
    plt.xlabel('Time [s]')
    plt.ylabel(r'Elevation [$^\circ$]')
    plt.xlim([0,max_time])
    plt.grid()
    
#linear fit
LF=np.polyfit(ppr,tau_nom,1) 

#%% Plots
plt.figure()
plt.errorbar(ppr,tau_nom,yerr=(tau_nom-tau_min,tau_max-tau_nom),fmt='.k',markersize=10,capsize=5)
plt.plot(ppr,LF[1]+LF[0]*ppr,'r')
plt.title('Source: ' +source)
plt.text(1000,1,r'$\Delta t='+str(np.round(LF[1],3))+'+'+str(np.round(LF[0]*10**4,3))+'\cdot 10^{-4}$ PPR')
plt.grid()
plt.xlabel('PPR')
plt.ylabel('$\Delta t$ [s]')
plt.ylim([0,1.25])
