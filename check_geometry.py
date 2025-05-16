# -*- coding: utf-8 -*-
'''
Quick visual check on geometry of scans
'''

import os
cd=os.path.dirname(__file__)
import sys
import yaml
import numpy as np
import linecache
import pandas as pd
from matplotlib import pyplot as plt

import matplotlib
plt.close('all')
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm' 
matplotlib.rcParams['font.size'] = 18

#%% Inputs
source_config='config.yaml'
source_data=os.path.join(cd,'data','awaken_multiPPT')#source of scans
source_scan=os.path.join(cd,'scans','250509_awaken_multiPPT')#source of scans
files={'User5_235_20250514_000016.hpl':'H04.w0.5.o20.xlsx',
       'User5_235_20250514_001000.hpl':'H05.w0.5.o20.xlsx',
       'User5_235_20250514_002000.hpl':'H06.w0.5.o20.xlsx'}
        
z_plot=500

#%% Initialization
with open(source_config, 'r') as fid:
    config = yaml.safe_load(fid)
    
#imports
sys.path.append(config['path_utils'])
import utils as utl



#%% Main
for file in files:
    f=os.path.join(source_data,file)
    
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
        
    scan_data=pd.read_excel(os.path.join(source_scan,files[file]))
    
    fig=plt.figure()
    ax=fig.add_subplot(111,projection='3d')
    for a,e in zip(azi,ele):
        rh=z_plot/utl.sind(e)
        plt.plot([0,rh*utl.cosd(e)*utl.cosd(90-a)],[0,rh*utl.cosd(e)*utl.sind(90-a)],[0,rh*utl.sind(e)],color='r',alpha=0.5)
        
    for a,e in zip(scan_data['Azimuth'].values,scan_data['Elevation'].values):
        rh=z_plot/utl.sind(e)
        plt.plot([0,rh*utl.cosd(e)*utl.cosd(90-a)],[0,rh*utl.cosd(e)*utl.sind(90-a)],[0,rh*utl.sind(e)],color='k',alpha=0.5)
    
    
