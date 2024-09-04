# -*- coding: utf-8 -*-
'''
Download data from Wind Data Hub
'''

import sys
import os
cd = os.path.dirname(__file__)

import warnings
import yaml

warnings.filterwarnings("ignore")

#%% Inputs
source_config=os.path.join(cd,'config.yaml')

#datset info
channel='awaken/sc1.lidar.z01.a0'
stime='20230306000000'#[%Y%m%d%H%M%S] start time UTC
etime='20230307000000'#[%Y%m%d%H%M%S] end time UTC
destination=''#absolute path to storage folder; if '' it will be ./data
file_format='nc'#format of files
extention='user2'#common attribute of files
MFA = False #multi factor authenticaiton needed? this is necessary for propietary datasets

#%% Initialization
with open(source_config, 'r') as fid:
    config = yaml.safe_load(fid)
    
#imports
sys.path.append(config['path_utils'])
sys.path.append(config['path_dap-py'])
import utils as utl
from doe_dap_dl import DAP

#%% Main
#connect to Wind Data Hub
a2e = DAP('a2e.energy.gov',confirm_downloads=False)
if MFA:
    a2e.setup_two_factor_auth(username=config['username'], password=config['password'])
else:
    a2e.setup_cert_auth(username=config['username'], password=config['password'])

# download data
if destination=='':
    destination = os.path.join(cd,'data',channel)
if os.path.exists(destination)==False:
    os.makedirs(destination)
_filter = {
    'Dataset': channel,
    'date_time': {
        'between': [stime,etime]
    },
    'file_type': 'nc',
    'ext1':extention
}
a2e.download_with_order(_filter, path=destination, replace=False)
