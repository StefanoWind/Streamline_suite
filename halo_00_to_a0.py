'''
Convert hpl Halo files to netCDF
'''
import os
cd=os.path.dirname(__file__)
import sys
sys.path.append('C:/Users/SLETIZIA/OneDrive - NREL/Desktop/Main/utils')
import utils as utl
import xarray as xr
import numpy as np
import yaml
import glob
import shutil
import re

#%% Inputs
source_config=os.path.join(cd,'config.yaml')
source='C:/Users/SLETIZIA/OneDrive - NREL/Desktop/Main/Technical/General/wdh-scripts/data/crosswind/nwtc.lidar.z01.00/*hpl'#source of hpl data files
site='nwtc'#site ID (sms=Siemens turbine)
instrument='lidar'
z_id='z01'#instrument index (if there are more instruments of the same type)
level='00'#input level of hpl (always 00)
save_path='C:/Users/SLETIZIA/OneDrive - NREL/Desktop/Main/Technical/General/wdh-scripts/data/crosswind'#source of hpl data files
replace=False
rename_files=False

#%% Functions

def rename(filename,site,instrument,z_id,level,save_path=None):
    
    if 'Stare' in filename:
        pattern = r"Stare_\d+_(\d{8})_(\d{2})(.*?)\.hpl"
        scan_type='stare'
    elif 'User' in filename:
        pattern = r"User\d{1}_\d+_(\d{8})_(\d{6})\.hpl"
        scan_type='user'
            
    match = re.search(pattern, os.path.basename(filename))
        
    date_part = match.group(1)  
    if len(match.group(2))<6:
        with open(filename, "r") as f:
            lines = []
            for line_num in range(11):
                lines.append(f.readline())
            metadata={}
            for line in lines:
                metaline = line.split(":")
                if "Start time" in metaline:
                    metadata["Start time"] = metaline[1:]
                else:
                    metadata[metaline[0]] = metaline[1]  # type: ignore
            time_part = match.group(2)+metadata["Start time"] [1]+metadata["Start time"] [2].split('.')[0]
    else:
        time_part = match.group(2)
        
    if save_path==None:
        directory=os.path.join('/'.join(os.path.dirname(filename).split('/')[:-1]),site+'.'+instrument+'.'+z_id+'.'+level)
        utl.mkdir(directory)
    else:
        directory=os.path.join(save_path,site+'.'+instrument+'.'+z_id+'.'+level)
    filename_out=site+'.'+instrument+'.'+z_id+'.'+level+'.'+date_part+'.'+time_part+'.'+scan_type+os.path.splitext(filename)[1]
    shutil.copyfile(filename, os.path.join(directory,filename_out))
    return os.path.join(directory,filename_out)

def read(filename):
    with open(filename, "r") as f:
        lines = []
        for line_num in range(11):
            lines.append(f.readline())

        # read metadata into strings
        metadata = {}
        for line in lines:
            metaline = line.split(":")
            if "Start time" in metaline:
                metadata["Start time"] = metaline[1:]
            else:
                metadata[metaline[0]] = metaline[1]  # type: ignore

        # convert some metadata
        num_gates = int(metadata["Number of gates"])  # type: ignore

        # Read some of the label lines
        for line_num in range(6):
            f.readline()

        # initialize arrays
        time = []
        azimuth = []
        elevation = []
        pitch = []
        roll = []
        doppler = []
        intensity = []
        beta = []

        while True:
            a = f.readline().split()
            if not len(a):  # is empty
                break

            time.append(float(a[0]))
            azimuth.append(float(a[1]))
            elevation.append(float(a[2]))
            pitch.append(float(a[3]))
            roll.append(float(a[4]))

            doppler.append([0] * num_gates)
            intensity.append([0] * num_gates)
            beta.append([0] * num_gates)

            for _ in range(num_gates):
                b = f.readline().split()
                range_gate = int(b[0])
                doppler[-1:][0][range_gate] = float(b[1])
                intensity[-1:][0][range_gate] = float(b[2])
                beta[-1:][0][range_gate] = float(b[3])

    # convert date to np.datetime64
    start_time_string = "{}-{}-{}T{}:{}:{}".format(
        metadata["Start time"][0][1:5],  # year
        metadata["Start time"][0][5:7],  # month
        metadata["Start time"][0][7:9],  # day
        "00",  # hour
        "00",  # minute
        "00.00",  # second
    )

    # find times where it wraps from 24 -> 0, add 24 to all indices after
    new_day_indices = np.where(np.diff(time) < -23)
    for new_day_index in new_day_indices[0]:
        time[new_day_index + 1 :] += 24.0

    start_time = np.datetime64(start_time_string)
    datetimes = [
        start_time + np.timedelta64(int(3600 * 1e6 * dtime), "us") for dtime in time
    ]

    dataset = xr.Dataset(
        {
            "Decimal time (hours)": (("time"), time),
            "Azimuth (degrees)": (("time"), azimuth),
            "Elevation (degrees)": (("time"), elevation),
            "Pitch (degrees)": (("time"), pitch),
            "Roll (degrees)": (("time"), roll),
            config['rws_name']: (("time", "range_gate"), doppler),
            "Intensity": (("time", "range_gate"), intensity),
            "Beta": (("time", "range_gate"), beta),
        },
        coords={"time": np.array(datetimes), "range_gate": np.arange(num_gates)},
        attrs={"Range gate length (m)": float(metadata["Range gate length (m)"])},  # type: ignore
    )

    # Save some attributes
    dataset.attrs[config['dr_name']] = float(
        metadata[config['dr_name']]  # type: ignore
    )
    dataset.attrs["Number of gates"] = float(metadata["Number of gates"])  # type: ignore
    dataset.attrs["Scan type"] = str(metadata["Scan type"]).strip()
    dataset.attrs[config['ppr_name']] = float(metadata["Pulses/ray"])  # type: ignore
    dataset.attrs["System ID"] = int(metadata["System ID"])  # type: ignore
    dataset.attrs["Filename"] = str(metadata["Filename"])[1:-5]
    dataset["distance"] = (
        "range_gate",
        dataset.coords["range_gate"].data * dataset.attrs["Range gate length (m)"]
        + dataset.attrs["Range gate length (m)"] / 2,
    )
    dataset["distance_overlapped"] = (
        "range_gate",
        dataset.coords["range_gate"].data * 1.5
        + dataset.attrs["Range gate length (m)"] / 2,
    )
    intensity = dataset.Intensity.data.copy()
    intensity[intensity <= 1] = np.nan
    dataset[config['snr_name']]=xr.DataArray(data=10 * np.log10(intensity - 1),coords={"time": np.array(datetimes), "range_gate": np.arange(num_gates)})


    # Dynamically add scan type and z-id (z02, z03, etc) to dataset metadata
    # loc_id, instrument, z02/z03, data level, date, time, scan type, extension
    raw_basename = filename.replace("\\", "/").rsplit("/")[-1]
    if ".z" in raw_basename:
        _, _, z_id, _, _, _, scan_type, _ = raw_basename.lower().split(".")
    else:  # local NREL tsdat-ing
        z_id = str(dataset.attrs["System ID"])
        scan_type = ""
        if "user" in raw_basename.lower():
            scan_type = raw_basename.split("_")[0].lower()
        elif "stare" in raw_basename.lower():
            scan_type = "stare"
        elif "vad" in raw_basename.lower():
            scan_type = "vad"
        elif "wind_profile" in raw_basename.lower():
            scan_type = "wind_profile"
        elif "rhi" in raw_basename.lower():
            scan_type = "rhi"

    valid_types = ["user", "stare", "vad", "wind_profile", "rhi"]
    if not any(valid_type in scan_type for valid_type in valid_types):
        raise NameError(f"Scan type '{scan_type}' not supported.")

    dataset.attrs["scan_type"] = scan_type
    dataset.attrs["z_id"] = z_id

    return dataset

#%% Initalization
files=glob.glob(source, recursive=True)

#config
with open(source_config, 'r') as fid:
    config = yaml.safe_load(fid)
    
#imports
sys.path.append(config['path_utils'])
import utils as utl

#%% Main
for f in files:
    # try:
    if rename_files:
        filename=rename(f,site,instrument,z_id,level,save_path)
    else:
        filename=f
    new_filename=os.path.join(os.path.dirname(filename.replace('.'+level,'.a0')),os.path.basename(filename).replace('.'+level+'.','.a0.').replace('hpl','nc'))
    if not os.path.exists(new_filename) or replace==True:
        dataset=read(filename)
        utl.mkdir(os.path.dirname(filename.replace('.'+level,'.a0')))
        dataset.to_netcdf(new_filename)
    # except:
    #     print(f+' failed')