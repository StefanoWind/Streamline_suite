# -*- coding: utf-8 -*-
"""
Simulator of Halo Streamline scanning head

Inputs:
    scan_file: name of scan text file (CSM or SSM)
    ppr: pulses per ray
    identifier: subtype of scan (e.g. overlapping, long range, short gate) as indicated in time info table
    mode: lidar model, as indicated in teme info table
    source_time: time info table file
"""
import pandas as pd
import time as tt
import re
import numpy as np
import sys
from matplotlib import pyplot as plt

def _round(x,resolution):
    return np.round(x/resolution)*resolution

def Halo_scan_sim(scan_file,ppr,identifier,model,source_time,dt=0.01,ang_tol=0.25,azi0=0,ele0=0):
    
    
    t0=tt.time()
    
    #read scan file
    with open(scan_file,'r') as fid:
        lines = fid.readlines()
    lines = [line.strip() for line in lines]
    
    if lines[0][:3]=='A.1':
        scan_mode='CSM'
    else:
        scan_mode='SSM'
    
    time_info=pd.read_excel(source_time,sheet_name=model)
    time_info=time_info.set_index('Scan settings')

    tp=time_info['Processing time'][scan_mode+ ' - '+identifier]
    ta=time_info['Acquisition time'][scan_mode+ ' - '+identifier]
    tau_s=tp+ta*ppr
    
    
    #zeroing
    P1=[]
    P2=[]
    S1=[]
    S2=[]
    A1=[]
    A2=[]
    
    time=[0]
    
    if scan_mode=='CSM':
        ppd1=500000/360#points per degree
        ppd2=250000/360#points per degree
        
        #extract kinematic parameters
        for l in lines:
            if len(l)>10:
                kin = [np.float64(v[1:]) for v in re.findall(r'=[-]?\d*\.?\d+', l)]
            
                A1=np.append(A1,kin[0])
                S1=np.append(S1,kin[1])
                P1=np.append(P1,kin[2])
                A2=np.append(A2,kin[3])
                S2=np.append(S2,kin[4])
                P2=np.append(P2,kin[5])
        
        #convert into S.I. units
        azi_range=np.append(azi0,-_round(P1/ppd1,ang_tol))
        ele_range=np.append(ele0,-_round(P2/ppd2,ang_tol))
        s1=S1*10/ppd1
        s2=S2*10/ppd2
        a1=A1*1000/ppd1
        a2=A2*1000/ppd2
        
        azi_all=azi_range[0]
        ele_all=ele_range[0]
    # elif scan_mode=='SSM':
    #     #inputs
    #     s1=time_info['Speed azimuth'][scan_mode+ ' - '+identifier]
    #     s2=time_info['Speed elevation'][scan_mode+ ' - '+identifier]
    #     a1=time_info['Acceleration azimuth'][scan_mode+ ' - '+identifier]
    #     a2=time_info['Acceleration elevation'][scan_mode+ ' - '+identifier]

        
    #scanning head movement        
    for i in range(len(P1)):
        azi_srt=azi_range[i]
        azi_end=azi_range[i+1]
        
        # if scan_mode=='SSM':#the azimuth in SSM fins the shortest path to the next location
        #     if th3-th0>180:
        #         th3=th3-360
        #     if th3-th0<-180:
        #         th3=th3+360
                
        ele_srt=ele_range[i]
        ele_end=ele_range[i+1]
        
        #zeroing
        t=[tau_s]
        dazi_dt=[0]
        azi=[azi_srt]
        t=[tau_s]
        dele_dt=[0]
        ele=[ele_srt]
        
        dazi_acc=np.min([s1[i]**2/(2*a1[i]+10**-10),np.abs(azi_end-azi_srt)/2])#azimuth span to accelrate/decelerate
        sign_azi=(azi_end-azi_srt)/np.abs(azi_end-azi_srt+10**-10)#azimuth direction
        
        dele_acc=np.min([s2[i]**2/(2*a2[i]+10**-10),np.abs(ele_end-ele_srt)/2])#elemuth span to accelrate/decelerate
        sign_ele=(ele_end-ele_srt)/np.abs(ele_end-ele_srt+10**-10)#elemuth direction
        
        for i_t in range(100000):
            t=np.append(t,t[-1]+dt)
            
            #check if it is time to decelerate (azimuth)
            if (azi_end-azi[-1])*sign_azi>dazi_acc:
                dazi_dt=np.append(dazi_dt,np.min([dazi_dt[-1]+a1[i]*dt,s1[i]]))#accelerate or keep max speed
            else:
                dazi_dt=np.append(dazi_dt,np.max([dazi_dt[-1]-a1[i]*dt,0]))#decelerate
            azi=np.append(azi,azi[-1]+dazi_dt[-1]*dt*sign_azi)
            
            #check if it is time to decelerate (elevation)
            if (ele_end-ele[-1])*sign_ele>dele_acc:
                dele_dt=np.append(dele_dt,np.min([dele_dt[-1]+a2[i]*dt,s2[i]]))#accelerate or keep max speed
            else:
                dele_dt=np.append(dele_dt,np.max([dele_dt[-1]-a2[i]*dt,0]))#decelerate
            ele=np.append(ele,ele[-1]+dele_dt[-1]*dt*sign_ele)

            #convergence
            if np.max([dazi_dt[-1],dele_dt[-1]])==0:
                break
        
        #round movement end
        azi[-1]=azi_end
        ele[-1]=ele_end
        
        # if scan_mode=='CSM':
        
        #append to global trajectory
        time=np.append(time,time[-1]+np.arange(t[0],t[-1],tau_s))
        azi_all=np.append(azi_all,np.interp(np.arange(t[0],t[-1],tau_s),t,azi))
        ele_all=np.append(ele_all,np.interp(np.arange(t[0],t[-1],tau_s),t,ele))
        
        #append end point
        time=np.append(time,time[-1]+tau_s)
        azi_all=np.append(azi_all,azi[-1])
        ele_all=np.append(ele_all,ele[-1])
        # elif scan_mode=='SSM':
        #     time=np.append(time,time[-1]+t[-1])
            
        ctr=i+1
        if np.floor(ctr/len(azi_range)*100)>np.floor((ctr-1)/len(azi_range)*100):
            est_time=(tt.time()-t0)/ctr*(len(azi_range)-ctr)
            sys.stdout.write('\r Scanning head simulator: '+str(np.floor(i/len(azi_range)*100).astype(int))+'% done, '+str(round(est_time))+' s left.') 
    sys.stdout.write('\r                                                                         ')
    sys.stdout.flush()
    
    return time,azi_all,ele_all
            
            
              
    