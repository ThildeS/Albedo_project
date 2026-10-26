

#%%
import os

import numpy as np
import pandas as pd
import astropy as astro
from datetime import datetime, timezone, timedelta
import time
from astral import moon
import pylunar
import math

import matplotlib
import matplotlib.pyplot as plt
print(matplotlib.get_backend())

# # # load functions from MT_functions
from Functions.JD2Lunar import *

# # check path
print('current working directory:', os.getcwd())

# # set data directory
data_dir = os.path.join(os.getcwd(), 'Data')

# # load txt file
filename = os.path.join(data_dir, 'BBSO_Albedo_anomaly.txt')
data = np.loadtxt(filename, skiprows=1)
print('data shape:', data.shape)

# # print first 5 rows of data
print('first 5 rows of data:')
print(data[:20])

# # plot the data inline in a Jupyter notebook
plt.figure(figsize=(12, 5))
plt.plot(data[:, 0], data[:, 1], label='Albedo Anomaly', color='cornflowerblue')
plt.axhline(y=0, color='black', lw=1)
plt.title('BBSO Albedo Anomaly Over Time')
plt.xlabel('Year')
plt.ylabel('Albedo Anomaly')
plt.legend()
plt.grid()
plt.show()

#%%

def decimal_year_to_utc(decimal_year):
    year = int(decimal_year)
    remainder = decimal_year - year
    
    # Check if the year is a leap year
    is_leap = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
    total_days = 366 if is_leap else 365
    
    # Calculate exact time delta from January 1st
    base_date = datetime(year, 1, 1, tzinfo=timezone.utc)
    utc_datetime = base_date + timedelta(days=remainder * total_days)
    return utc_datetime

def get_julian_datetime(date_to_convert):
    utc_datetime = decimal_year_to_utc(date_to_convert)
    year = utc_datetime.year
    month = utc_datetime.month
    day = utc_datetime.day
    hour = utc_datetime.hour
    minute = utc_datetime.minute
    second = utc_datetime.second

    date = datetime(year, month, day, hour, minute, second)

    # Ensure correct format
    if not isinstance(date, datetime):
        raise TypeError('Invalid type for parameter "date" - expecting datetime')
    elif date.year < 1801 or date.year > 2099:
        raise ValueError('Datetime must be between year 1801 and 2099')

    # Perform the calculation
    julian_datetime = 367 * date.year - int((7 * (date.year + int((date.month + 9) / 12.0))) / 4.0) + int(
        (275 * date.month) / 9.0) + date.day + 1721013.5 + (
                          date.hour + date.minute / 60.0 + date.second / math.pow(60,
                                                                                  2)) / 24.0 - 0.5 * math.copysign(
        1, 100 * date.year + date.month - 190002.5) + 0.5

    return julian_datetime

def moon_phase(decimal_year):
    utc_datetime = decimal_year_to_utc(decimal_year)
    year = utc_datetime.year
    month = utc_datetime.month
    day = utc_datetime.day
    
    # Calculate moon phase using pylunar
    moon_phase_value = moon.phase(datetime.date(year, month, day))
    
    return moon_phase_value

def moon_phase_name(decimal_year):
    utc_datetime = decimal_year_to_utc(decimal_year)
    year = utc_datetime.year
    month = utc_datetime.month
    day = utc_datetime.day
    
    # Calculate moon phase name using pylunar
    mi = pylunar.MoonInfo((34, 15, 30), (116, 55, 17))  # BBSO coordinates (latitude, longitude)
    mi.update((year, month, day))
    phase_name = mi.phase_name()
    
    return phase_name


# turn data into pandas dataframe
df = pd.DataFrame(data, columns=['Decimal_Year', 'Albedo_Anomaly'])

# add columns for UTC datetime, moon phase, and moon phase name
df.insert(1, 'UTC_datetime', df['Decimal_Year'].apply(decimal_year_to_utc))

# add julian dates
df.insert(2, 'Julian_date', df['Decimal_Year'].apply(get_julian_datetime))

print('first 5 rows of dataframe:')
print(df.head())

#df.insert(3, ['phase_angle', 'sublat', 'sublon', 'illum'], df['Julian_date'].apply(lunar_phase_and_sublunar_from_jd))

values = lunar_phase_and_sublunar_from_jd(df['Julian_date'][0])
print(values)

#new_cols = ['phase_angle', 'sublat', 'sublon', 'illum']
new_cols = ['phase_angle', 'signed_phase', 'sublat', 'sublon', 'illum']
values_to_append = df['Julian_date'].apply(lunar_phase_and_sublunar_from_jd)
temp_df = pd.DataFrame(values_to_append.tolist(), columns=new_cols)
df = pd.concat([df, temp_df], axis=1)

# save dataframe to txt file
# output_filename = os.path.join(data_dir, 'BBSO_Data.txt')
# df.to_csv(output_filename, sep='\t', index=False)

# remove data from before 13 June 2015
print('\nRemove data')
BBSO_EPIC = df[df['UTC_datetime'] > '2015-06-13']
print(BBSO_EPIC.head())
print(f'Number of days : {len(BBSO_EPIC)}')

plt.figure(figsize=(12, 5))
plt.plot(BBSO_EPIC['Decimal_Year'], BBSO_EPIC['Albedo_Anomaly'], label='Albedo Anomaly', color='cornflowerblue')
plt.axhline(y=0, color='black', lw=2)
plt.title('BBSO Albedo Anomaly EPIC data overlap')
plt.xlabel('Year')
plt.ylabel('Albedo Anomaly')
plt.legend()
plt.grid()
plt.show()

# output_filename = os.path.join(data_dir, 'BBSO_EPIC_data.txt')
# BBSO_EPIC.to_csv(output_filename, sep='\t', index=False)

# make list of days used taking the UTC format at writing is epic_1b__YYYYMMDD* where YYYY is the year, MM is the month, and DD is the day
days_used = BBSO_EPIC['UTC_datetime'].dt.strftime('epic_1b_%Y%m%d*').tolist()
print(days_used[:10])

# remove the string '' markers from the list
days_used = [day.replace("'", "") for day in days_used]
print(days_used[:10])

# save list of days used to txt file
# output_filename = os.path.join(data_dir, 'BBSO_EPIC_days_used.txt')
# with open(output_filename, 'w') as f:
#     for day in days_used:
#         f.write(f"{day}\n")



#%%

import h5py

#Open the H5 file in read mode
file = data_dir + '/' + 'epic_1b_20150621163526_04.h5'

f = h5py.File(file, 'r')
list(f.keys())
print(list(f.keys()))
band_name = 'Band680nm'
B = f[band_name]

print(B)
print(B.keys())

gloc = B['Geolocation']
img = B['Image']
pt = B['PixelType']
print(img)
print(img.dtype)

# get image information
d_img = dict(img.attrs)

print('\nPrinting image file attributes \n')
for keys, values in d_img.items():
    print(f'{keys}: {values}')

plt.figure()
plt.imshow(img[:].T)
cbar = plt.colorbar()
plt.title(band_name)
cbar.set_label('counts/second')
plt.show()

print(d_img.get('time'))
print(f['Band443nm']['Image'].attrs['time'])


time_obj = f['Band443nm']['Image'].attrs['time']
time_str = time_obj.decode('utf-8')

# 2. Parse into a datetime object
dt = datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S')

# 3. Separate into date and time objects
date_val = dt.date()   # Output: datetime.date(2015, 6, 21)
clock_val = dt.time()  # Output: datetime.time(16, 30, 38)

# If you just want them as separate strings instead:
date_str = dt.strftime('%Y-%m-%d') # '2015-06-21'
clock_str = dt.strftime('%H:%M:%S') # '16:30:38'
print(clock_str)
#date = 

# d_gloc = dict(gloc.attrs)
# print(d_gloc)

# open geolocation group
print(gloc)

# plot the geolocation data



# %%
