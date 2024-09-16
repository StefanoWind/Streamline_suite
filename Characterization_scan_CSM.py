# -*- coding: utf-8 -*-
'''
Design scan to characterize sampling time in CSM mode
'''

import os
cd=os.path.dirname(__file__)
import sys
import yaml
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from datetime import datetime 

plt.close('all')
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm' 
matplotlib.rcParams['font.size'] = 14

#%% Inputs
source_config='config.yaml'
narrow_scan_amplitude=40#[deg] amlpitude of narrow scans
N_speeds=5#number of tested speeds

max_S=5000#[10's of pulses per second] motor maximum speed
max_A=50#[1000's of pulses per second^2] motor maximum acceleration
ppd1=500000/360#[pulses per degree] firing rate in azimuth
ppd2=250000/360#[pulses per degree] firing rate elevation

#%% Initialization
#config
with open(source_config, 'r') as fid:
    config = yaml.safe_load(fid)
    
#imports
sys.path.append(config['path_utils'])
import utils as utl

#define scanning speeds
speeds=utl.round(np.linspace(max_S/10,max_S,N_speeds),100)

#zeroing
P1=[]
S1=[]
A1=[]
P2=[]
S2=[]
A2=[]

L=''

#%% Main

#narrow PPIs
for S in speeds:
    #go to start
    P1=np.append(P1,0)
    S1=np.append(S1,max_S)
    A1=np.append(A1,max_A)
    P2=np.append(P2,0)
    S2=np.append(S2,max_S)
    A2=np.append(A2,max_A)
    
    #move at test speed
    P1=np.append(P1,-narrow_scan_amplitude*ppd1)
    S1=np.append(S1,S)
    A1=np.append(A1,max_A)
    P2=np.append(P2,0)
    S2=np.append(S2,0)
    A2=np.append(A2,0)
    
#narrow RHIs
for S in speeds:
    #go to start
    P1=np.append(P1,0)
    S1=np.append(S1,max_S)
    A1=np.append(A1,max_A)
    P2=np.append(P2,0)
    S2=np.append(S2,max_S)
    A2=np.append(A2,max_A)
    
    #move at test speed
    P1=np.append(P1,0)
    S1=np.append(S1,0)
    A1=np.append(A1,0)
    P2=np.append(P2,-narrow_scan_amplitude*ppd2)
    S2=np.append(S2,S)
    A2=np.append(A2,max_A)
    
#wide PPIs
for S in speeds:
    #go to start
    P1=np.append(P1,0)
    S1=np.append(S1,max_S)
    A1=np.append(A1,max_A)
    P2=np.append(P2,0)
    S2=np.append(S2,max_S)
    A2=np.append(A2,max_A)
    
    #move at test speed
    P1=np.append(P1,-359*ppd1)
    S1=np.append(S1,S)
    A1=np.append(A1,max_A)
    P2=np.append(P2,0)
    S2=np.append(S2,0)
    A2=np.append(A2,0)
    

#wide RHIs
for S in speeds:
    #go to start
    P1=np.append(P1,0)
    S1=np.append(S1,max_S)
    A1=np.append(A1,max_A)
    P2=np.append(P2,0)
    S2=np.append(S2,max_S)
    A2=np.append(A2,max_A)
    
    #move at test speed
    P1=np.append(P1,0)
    S1=np.append(S1,0)
    A1=np.append(A1,0)
    P2=np.append(P2,-179*ppd2)
    S2=np.append(S2,S)
    A2=np.append(A2,max_A)
    
#write scan file
for p1,p2,s1,s2,a1,a2 in zip(P1,P2,S1,S2,A1,A2):
    
    l='A.1=%.0f'%a1+',S.1=%.0f'%s1+',P.1=%.0f'%p1+'*A.2=%.0f'%a2+',S.2=%.0f'%s2+',P.2=%.0f'%p2+'\nW=0\n'
    L+=l

#%% Output
utl.mkdir('scans')
name=os.path.join(cd,'scans',datetime.strftime(datetime.now(),'%y%m%d.%H%M%S')+'.CSM_test.txt')

with open(name,'w') as fid:
    fid.write(L)
    fid.close()