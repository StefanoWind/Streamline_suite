# -*- coding: utf-8 -*-
'''
Estimate noise from stare scan using the autocorrelatin method (Lenschow et el, 2000)
'''

import os
cd=os.path.dirname(__file__)
import sys
import numpy as np
import xarray as xr
import yaml
from matplotlib import pyplot as plt
import glob
from scipy import stats
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

import matplotlib
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm' 
matplotlib.rcParams['font.size'] = 14
plt.close('all')

#%% Inputs
source_config=os.path.join(cd,'config.yaml')

if len(sys.argv)==1:
    source=os.path.join(cd,'data/test/*202411*nc')
    output_name='nwtc.lidar.z01'#name of output files
else:
    source=sys.argv[1]
    output_name=sys.argv[2]

#noise estimation
N_lags=100#number of lags of ACF
max_nsi=0.5#maximum ratio of standard deviation standard deviation to mean standard deviation
DT=600#[m] averaging time
min_time_bins=5#number of time bins to check non-stationarity
p_value=0.05#p_value for bootstrap
bins_snr=np.arange(-40.5,10.6)#[dB] bins in snr

#qc
rmin=96#[m] blind zone of the lidar
rmax=3000#[m] max range 
rws_max=40#[m/s] maximum rws discarded prior to the ACF estimation
min_noise=10**-10 #[m/s] minimum noise level
min_Nt=100#minimum number of samples
prec_rws_thresh=-1#[m/s] rws threshold to flag precipitation (smaller than)
prec_snr_thresh=-10#[dB] snr threshold to flag precipitation (higher than)
prec_r_thresh=1000#[m] range threshold to flag precipitation (higher than)

#%% Initialization

#config
with open(source_config, 'r') as fid:
    config = yaml.safe_load(fid)
    
#imports
sys.path.append(config['path_utils'])
import utils as utl

#files
files=sorted(glob.glob(os.path.join(cd, source)))

#stats
snr_avg=utl.mid(bins_snr)

#zeroing
ACF_all={}
snr_all={}
noise={}
noise_avg={}
noise_top={}
noise_low={}

#create folders
os.makedirs(os.path.join(cd,'data',output_name),exist_ok=True)
os.makedirs(os.path.join(cd,'figures',output_name),exist_ok=True)

#%% Main
for f in files:
    print(f'Processing {os.path.basename(f)}', flush=True)
    with xr.open_dataset(f) as Data:
        
        #sort time
        Data=Data.sortby('time')
        
        #get rid of duplicate timestamps
        _, index = np.unique(Data["time"], return_index=True)
        Data = Data.isel(time=index)
        
        #select physical distance
        if 'overlapping' in Data.attrs['Scan type']:
            r0=Data['distance_overlapped'].values
        else:
            r0=Data['distance'].values
        sel_r=(r0>rmin)*(r0<rmax)
        r=r0[sel_r]
        
        #remove and interpolate large rws values
        tnum=np.array(Data['time'].values.tolist())/10**9
        Data[config['rws_name']]=Data[config['rws_name']].where(np.abs(Data[config['rws_name']])<rws_max).interpolate_na(dim='time')
        
        #interpolate time on uniform linear space
        PF=np.polyfit(np.arange(len(tnum)), tnum, 1)
        tnum_uni=PF[1]+PF[0]*np.arange(len(tnum))
        time_uni=tnum_uni*10**9*np.timedelta64(1,'ns')+np.datetime64('1970-01-01T00:00:00')
        Data_uni=Data.interp(time=time_uni)
        
        #extract variables
        rws=np.array(Data_uni[config['rws_name']])[:,sel_r].T
        snr=np.array(Data_uni[config['snr_name']])[:,sel_r].T

        ppr=int(Data.attrs[config['ppr_name']])
        dr=int(Data.attrs[config['dr_name']])
    
    #precipitation check
    prec=np.sum((np.nanmean(rws[r>prec_r_thresh,:],axis=1)<prec_rws_thresh)*(np.nanmean(snr[r>prec_r_thresh,:],axis=1)>prec_snr_thresh))
    
    if prec>0:
        print(f'Possible precipitation detected, skipping {os.path.basename(f)}', flush=True)
        continue
    
    #stationarity check (based on stddev of binned stdev)
    bin_tnum_nsi=np.arange(tnum[0],tnum[-1]+DT/2,DT)
    if len(bin_tnum_nsi)<min_time_bins+1:
        bin_tnum_nsi=np.linspace(tnum[0],tnum[-1],min_time_bins+1)
    rws_std=[]
    snr_std=[]
    for t1,t2 in zip(bin_tnum_nsi[:-1],bin_tnum_nsi[1:]):
        sel_t=(tnum>=t1)*(tnum<=t2)
        rws_std=utl.vstack(rws_std,np.nanstd(rws[:,sel_t],axis=1))
        snr_std=utl.vstack(snr_std,np.nanstd(snr[:,sel_t],axis=1))
        
    rws_nsi=np.nanstd(rws_std,axis=0)/np.nanmean(rws_std,axis=0)#non-stationarity index for rws
    snr_nsi=np.nanstd(snr_std,axis=0)/np.nanmean(snr_std,axis=0)#non-stationarity index for snr
    
    rws_qc=rws.copy()
    rws_qc[rws_nsi>max_nsi,:]=np.nan
    rws_qc[snr_nsi>max_nsi,:]=np.nan
    
    setup='ppr='+str(ppr)+'.dr='+str(dr)
    
    #noise estimation
    bin_tnum=np.arange(tnum[0],tnum[-1]+DT/2,DT)
    for t1,t2 in zip(bin_tnum[:-1],bin_tnum[1:]):
        sel_t=(tnum>=t1)*(tnum<=t2)
        Nt=np.sum(sel_t)
        if Nt>min_Nt:
            ACF=np.zeros((len(r),N_lags))
            rws_avg=np.tile(np.nanmean(rws_qc[:,sel_t],axis=1),(Nt,1)).T
            rws_det=rws_qc[:,sel_t]-rws_avg
            for i_r in range(len(r)):
                conv=np.correlate(rws_det[i_r,:],rws_det[i_r,:], mode='full')
                N=np.correlate(np.zeros(Nt)+1, np.zeros(Nt)+1, mode='full')
                ACF[i_r,:]=conv[Nt-1:Nt-1+N_lags]/N[Nt-1:Nt-1+N_lags]
                
                #check on variance vs 0-lag ACF
                if np.abs(ACF[i_r,0]-np.nanvar(rws_det[i_r,:]))>10**-10:
                    raise ValueError('Variance mismatch')
            
            #initialize structures for new lidar setup
            if setup not in ACF_all:
                ACF_all[setup]=[]
                snr_all[setup]=[]
    
            ACF_all[setup]=utl.vstack(ACF_all[setup],ACF)
            snr_all[setup]=np.append(snr_all[setup],np.nanmean(snr[:,sel_t],axis=1))

#overall statistics
for s in ACF_all.keys():
    #linear extrapolation of noise
    ACF_id=ACF_all[s].copy()
    ACF_id[:,0]=(ACF_all[s][:,1]+(ACF_all[s][:,1]-ACF_all[s][:,2]))
    noise[s]=(ACF_all[s][:,0]-ACF_id[:,0])**0.5
    noise[s][ACF_all[s][:,1]<ACF_all[s][:,2]]=np.nan
    noise[s][noise[s]<min_noise]=np.nan
    
    #statistics
    noise_avg[s]=10**stats.binned_statistic(snr_all[s],np.log10(noise[s]),statistic=lambda x:utl.filt_mean(x),                       bins=bins_snr)[0]
    noise_low[s]=10**stats.binned_statistic(snr_all[s],np.log10(noise[s]),statistic=lambda x:utl.filt_BS_mean(x,p_value/2*100),      bins=bins_snr)[0]
    noise_top[s]=10**stats.binned_statistic(snr_all[s],np.log10(noise[s]),statistic=lambda x:utl.filt_BS_mean(x,(1-p_value/2)*100),  bins=bins_snr)[0]

    snr_hist=np.histogram(snr_all[s][noise[s]>0],bins=bins_snr)[0]
    
    #Output
    Output=pd.DataFrame()
    Output['SNR [dB]']=snr_avg
    Output['Noise StDev [m/s]']=noise_avg[s]
    Output[f'Noise StDev ({p_value/2*100}% percentile) [m/s]']=noise_low[s]
    Output[f'Noise StDev ({(1-p_value/2)*100}% percentile) [m/s]']=noise_top[s]
    Output['Occurrence']=snr_hist
    
    with pd.ExcelWriter(os.path.join(cd,'data',output_name,output_name+s+'.snr.noise.cutoff'+str(rws_max)+'.xlsx')) as writer:
        Output.to_excel(writer, sheet_name=s, index=False)
    
    #Plots
    fig=plt.figure(figsize=(18,9))
    main_ax = fig.add_axes([0.1, 0.3, 0.6, 0.6]) 
    plt.semilogy(snr_all[s],noise[s],'.k',alpha=0.05, markersize=10)
    plt.errorbar(snr_avg,noise_avg[s],[noise_avg[s]-noise_low[s],noise_top[s]-noise_avg[s]],color='r')
    reals=~np.isnan(noise_avg[s])
    plt.plot(snr_avg[reals],noise_avg[s][reals],'.-r',markersize=15)
    plt.grid()
    plt.title('Noise estimation for '+s+' based on '+str(len(ACF_all[s][:,0]))+' stares')
    plt.xlabel(r'SNR [dB]')
    plt.ylabel('Measured noise Standard Deviation [m s$^{-1}$]')
    plt.xlim([-30,-5])
    plt.xticks(np.arange(-30,-9,5))
    plt.ylim([0.01,30])
    
    utl.mkdir(os.path.join(cd,'figures'))
    plt.savefig(os.path.join(cd,'figures',output_name,output_name+s+'.snr.noise.cutoff'+str(rws_max)+'.png'))
    plt.close()

    

