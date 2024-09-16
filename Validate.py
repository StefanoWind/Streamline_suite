# -*- coding: utf-8 -*-
'''
Validate scanning head simulator
'''

import os
cd=os.path.dirname(__file__)
import sys
import yaml
import numpy as np
import linecache
import glob 
from matplotlib import pyplot as plt
import Halo_scan_sim as HSS
import matplotlib
plt.close('all')
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm' 
matplotlib.rcParams['font.size'] = 16

#%% Inputs
source_config='config.yaml'
source=os.path.join(cd,'data','CSM_test','*.hpl')#source of scans
scan_file=os.path.join(cd,'scans/240828.114620.CSM_test.txt')#scan file
source_time=os.path.join(cd,'data/Halo_time_info.xlsx')#time information sheet
identifier='no-overlapping'
model='XR (Crosswind)'# lidar model
ang_tol=0.25#[deg] angular tolerance

#graphics
max_time=410#[s] maixmum time for plot

#%% Initialization
with open(source_config, 'r') as fid:
    config = yaml.safe_load(fid)
    
#imports
sys.path.append(config['path_utils'])
import utils as utl

#file names
files=glob.glob(os.path.join(cd,'data',source))

#%% Main
for f in files:
    #zeroing
    time=[]
    ele=[]
    azi=[]
    
    #read real data
    N_tot = sum(1 for line in open(f))
    Nr=int(linecache.getline(f, 3).split(':')[1])
    ppr=int(linecache.getline(f, 6).split(':')[1])
    
    line_number=18
    ctr=0
    while line_number<N_tot:
        l=linecache.getline(f, line_number).split(' ')
        l_strip=[]
        for ll in l:
            if len(ll)>0:
                l_strip.append(ll)
        
        time=np.append(time,np.float64(l_strip[0])*3600) 
        azi=np.append(azi,np.float64(l_strip[1]))
        ele=np.append(ele,np.float64(l_strip[2]))          
        ctr+=1
        line_number=18+ctr*(Nr+1)         
    
    #fix coordinates
    time=time-time[0]
    azi=utl.round(azi,ang_tol)%360
    ele=utl.round(ele,ang_tol)%360
        
    #run scanning head simulator
    time_sim,azi_sim,ele_sim=HSS.Halo_scan_sim(scan_file,ppr,identifier,model,source_time)

    # plots
    plt.figure(figsize=(18,8))
    plt.subplot(1,2,1)
    plt.plot(time,azi,'.-k',label='Data')
    plt.plot(time_sim,azi_sim,'o-r',fillstyle='none',markersize=3,label='Model')
    plt.grid()
    plt.xlim([0,max_time])
    plt.xlabel('Time [s]')
    plt.ylabel(r'$\theta$ [$^\circ$]')
    plt.title(' '*120+'Real trajectory from '+os.path.basename(f)+', \n'+' '*120+' simulated trajectory based on '+os.path.basename(scan_file)+', PPR='+str(ppr))
        
    plt.subplot(1,2,2)
    plt.plot(time,ele,'.-k',label='Data')
    plt.plot(time_sim,ele_sim,'o-r',fillstyle='none',markersize=3,label='Model')
    plt.grid()
    plt.xlim([0,max_time])
    plt.xlabel('Time [s]')
    plt.ylabel(r'$\beta$ [$^\circ$]')
    plt.legend()
    plt.tight_layout()
    
    
    
    
