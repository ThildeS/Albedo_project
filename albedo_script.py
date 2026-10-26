
#%%
import os
import numpy as np
import pandas as pd
import astropy as astro
from datetime import datetime, timezone, timedelta, time
#import time
import math
import h5py
h5py.enable_ipython_completer()

import matplotlib
import matplotlib.pyplot as plt
print(matplotlib.get_backend())

# # # load functions from functions
#from Functions.JD2Lunar import *

### start of the code

# # check path
print('current working directory:', os.getcwd())

# # set data directory
data_dir = os.path.join(os.getcwd(), 'Data')
image_dir = os.path.join(os.getcwd(), 'Data')

# # load txt file
filename = os.path.join(data_dir, 'BBSO_EPIC_data.txt')
column_names = ['Decimal_Year', 'UTC_datetime', 'Julian_date', 'Albedo_Anomaly', 'phase_angle', 'signed_phase', 'sublat', 'sublon', 'illum']
data = pd.read_csv(filename, sep = "\t", names = column_names, skiprows=1)
print('data shape:', data.shape)
print(data.head())

### functions ###


def Get_image(f, band):
    image = f[band]['Image']
    image = np.fliplr(np.rot90(image, -1))

    # mask out background
    mask_grid = f[band]['Geolocation']['Earth']['Mask']
    #print('Mask grid shape:', mask_grid.shape)
    mask_grid = np.fliplr(np.rot90(mask_grid, -1))
    mask = mask_grid > 0
    image_masked = np.ma.masked_array(image, mask=~mask)

    # mask the reflectance image where sza > 80 degrees
    sza = f[band]['Geolocation']['Earth']['SunAngleZenith']
    sza = np.fliplr(np.rot90(sza, -1))
    sza_mask = (sza > 78)
    image_masked = np.ma.masked_array(image_masked, mask=sza_mask)

    return image_masked

def Get_grid(f, band):
    lat_grid = f[band]['Geolocation']['Earth']['Latitude']
    lat_grid = np.fliplr(np.rot90(lat_grid, -1))

    lon_grid = f[band]['Geolocation']['Earth']['Longitude']
    lon_grid = np.fliplr(np.rot90(lon_grid, -1))
    return lat_grid, lon_grid

def Area_pixels(f, band, region_tag):
    lat_grid, lon_grid = Get_grid(f, band)
    if region_tag == 'Atlantic':
        # mask the Atlantic region
        region_mask = (lat_grid >= -69) & (lat_grid <= 67) & (lon_grid >= -85) & (lon_grid <= 60)
    else:
        region_mask = (
        (lat_grid >= -69) & (lat_grid <= 67) & 
        (((lon_grid > -180) & (lon_grid <= -150)) | ((lon_grid >= 70) & (lon_grid <= 180)))
        )
    earth_mask = f[band]['Geolocation']['Earth']['Mask']
    earth_pixels = np.nansum(earth_mask)
    region_pixels = np.nansum(region_mask)
    
    return region_pixels, earth_pixels

def Area_mask(f, band, region_tag):
    image = Get_image(f, band)
    lat_grid, lon_grid = Get_grid(f, band)
    if region_tag == 'Atlantic':
        # mask the Atlantic region
        region_mask = (lat_grid >= -69) & (lat_grid <= 67) & (lon_grid >= -85) & (lon_grid <= 60)
    else:
        region_mask = (
        (lat_grid >= -69) & (lat_grid <= 67) & 
        (((lon_grid > -180) & (lon_grid <= -150)) | ((lon_grid >= 70) & (lon_grid <= 180)))
        )
    image_masked = np.ma.masked_array(image, mask=~region_mask)
    return image_masked

def get_info(f):
    # centroid coordinates
    cent_lat = float(np.asarray(f['Band443nm']['Geolocation']['Earth'].attrs['centroid_center_latitude']).squeeze())
    cent_lon = float(np.asarray(f['Band443nm']['Geolocation']['Earth'].attrs['centroid_center_longitude']).squeeze())

    # time
    UTC_time = f['Band680nm']['Image'].attrs['time']
    UTC_time = UTC_time.decode('utf-8')
    dt = datetime.strptime(UTC_time, '%Y-%m-%d %H:%M:%S')
    date_str = dt.strftime('%Y-%m-%d')
    clock_str = dt.strftime('%H:%M:%S')

    return cent_lat, cent_lon, date_str, clock_str


def Plot():
    print(f"Plotting: {file}")
    # Create one figure with three panels
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, band in zip(axes, bands):
        image = Area_mask(f, band, region_tag)
        image_ref = image * RCF[band]

        im = ax.imshow(image_ref)
        ax.set_title(band)
        fig.colorbar(im, ax=ax)

    fig.suptitle(f"{file}\n{image_date} {clock_str} | {region_tag}")
    fig.tight_layout()
    plt.show()
    plt.close(fig)

    print(f"Finished plotting: {file}")


# %%
albedo_data = []
column_names = [
    "Filename",
    "Date",
    "Time",
    "Centroid_lat",
    "Centroid_lon",
    "Region",
    "Region_pixels",
    "Coverage",
    "Albedo_443nm",
    "Albedo_551nm",
    "Albedo_680nm",
]
df = pd.DataFrame(albedo_data, columns=column_names)

bands = ['Band443nm', 'Band551nm', 'Band680nm']
RCF = {'Band443nm': 8.340e-06,
        'Band551nm': 6.660e-06,
        'Band680nm': 9.300e-06}

########### loop through the data  #################

for index, row in data.iterrows():
    # get date string
    date = row['UTC_datetime'].split(' ')[0]
    date_str = date.replace('-', '')

    # get the lunar phase
    signed_phase = row['signed_phase']
    region_tag = ''

    if signed_phase > 0: # the waning phase = Atlantic region
        region_tag = 'Atlantic'
    else: # the waxing phase = Astralasian region
        region_tag = 'Australasian'

    # read only images from that date within the time window using the file name epic_1b_20160211132212_04.h5
    for file in sorted(os.listdir(image_dir)):
        if file.startswith(f'epic_1b_{date_str}') and file.endswith('.h5'):
            file_date_str = file.split('_')[2]
            datetime_obj = datetime.strptime(file_date_str, '%Y%m%d%H%M%S')
            time_obj = datetime_obj.time() # Extract time as a datetime.time object
            print(f"\nChecking {file}")

            if region_tag == 'Atlantic':
                if time(10, 0, 0) <= time_obj <= time(15, 0, 0):
                    print(f"SELECTED: {file} | time: {time_obj}, region: {region_tag}")
                    print(f"START processing")
                    try:
                        with h5py.File(os.path.join(image_dir, file), 'r') as f:
                            cent_lat, cent_lon, image_date, clock_str = get_info(f)
                            band_albedo = {}
                            for band in bands:
                                print(f"Processing {band}")
                                image = Area_mask(f, band, region_tag)
                                image_ref = image * RCF[band]

                                # calculate albedo
                                valid_pixels = image_ref.compressed()
                                if valid_pixels.size > 0:
                                    band_albedo[band] = np.nanmean(valid_pixels)
                                else:
                                    band_albedo[band] = np.nan
                                print(f"{band} completed")

                            # Find the regional coverage
                            region_pixels, earth_pixels = Area_pixels(f, "Band443nm", region_tag)
                            coverage = region_pixels / earth_pixels

                        print(f"FINISHED processing: {file}, Saving results")
                        new_row = [
                            file,
                            image_date,
                            clock_str,
                            cent_lat,
                            cent_lon,
                            region_tag,
                            region_pixels,
                            coverage,
                            band_albedo["Band443nm"],
                            band_albedo["Band551nm"],
                            band_albedo["Band680nm"],
                        ]
                        albedo_data.append(new_row)
                        print(new_row)

                    except Exception as e:
                        print(f"ERROR in {file}: {type(e).__name__}: {e}")
                        continue
                else:
                    print(f"Outside time window")

            else:
                if time(0, 0, 0) <= time_obj <= time(5, 0, 0):
                    print(f"SELECTED: {file} | time: {time_obj}, region: {region_tag}")
                    print(f"START processing")
                    try:
                        with h5py.File(os.path.join(image_dir, file), 'r') as f:
                            cent_lat, cent_lon, image_date, clock_str = get_info(f)
                            band_albedo = {}
                            for band in bands:
                                print(f"Processing {band}")
                                image = Area_mask(f, band, region_tag)
                                image_ref = image * RCF[band]

                                # calculate albedo
                                valid_pixels = image_ref.compressed()
                                if valid_pixels.size > 0:
                                    band_albedo[band] = np.nanmean(valid_pixels)
                                else:
                                    band_albedo[band] = np.nan
                                print(f"{band} completed")

                            # Find the regional coverage
                            region_pixels, earth_pixels = Area_pixels(f, "Band443nm", region_tag)
                            coverage = region_pixels / earth_pixels

                        print(f"FINISHED processing: {file}, Saving results")
                        new_row = [
                            file,
                            image_date,
                            clock_str,
                            cent_lat,
                            cent_lon,
                            region_tag,
                            region_pixels,
                            coverage,
                            band_albedo["Band443nm"],
                            band_albedo["Band551nm"],
                            band_albedo["Band680nm"],
                        ]
                        albedo_data.append(new_row)
                        print(new_row)

                    except Exception as e:
                        print(f"ERROR in {file}: {type(e).__name__}: {e}")
                        continue
                else:
                    print(f"Outside time window")

df = pd.DataFrame(albedo_data, columns=column_names)
print(df.head())

output_filename = os.path.join(data_dir, 'EPIC_albedo.txt')
df.to_csv(output_filename, sep='\t', index=False)


# %%
