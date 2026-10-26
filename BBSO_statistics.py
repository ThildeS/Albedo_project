#%% BBSO statistics

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


# # check path
print('current working directory:', os.getcwd())

# # set data directory
data_dir = os.path.join(os.getcwd(), 'Data')

# # load txt file
filename = os.path.join(data_dir, 'BBSO_Data.txt')
column_names = ['Decimal_Year', 'UTC_datetime', 'Julian_date', 'Albedo_Anomaly', 'phase_angle', 'signed_phase', 'sublat', 'sublon', 'illum']
data = pd.read_csv(filename, sep = "\t", names = column_names, skiprows=1)
print('data shape:', data.shape)

print(data.head())

# seperate the data into two parts: illum < 0.3 and illum >= 0.3
data_low = data[data['illum'] < 0.3]
data_high = data[data['illum'] >= 0.3]

# %% Plot Albedo anomaly vs illumination fraction

plt.figure(figsize=(12, 5))
plt.scatter(data_low['illum'], data_low['Albedo_Anomaly'], label='illum < 0.3', color='cornflowerblue', alpha=0.5)
plt.scatter(data_high['illum'], data_high['Albedo_Anomaly'], label='illum >= 0.3', color='orange', alpha=0.5)
plt.axhline(y=0, color='black', lw=1)
plt.xlabel('Illumination Fraction')
plt.ylabel('Albedo Anomaly')
plt.title('BBSO Albedo Anomaly vs Illumination Fraction for 1998-2017')
plt.legend()
plt.show()

# make boxplot of Albedo anomaly vs illumination fraction
plt.figure(figsize=(12, 5))
plt.boxplot([data_low['Albedo_Anomaly'], data_high['Albedo_Anomaly']], labels=['illum < 0.3', 'illum >= 0.3'])
plt.axhline(y=0, color='black', lw=1)
plt.ylabel('Albedo Anomaly')
plt.title('BBSO Albedo Anomaly vs Illumination Fraction for 1998-2017')
plt.show()

# plot illumination fraction vs decimal year
plt.figure(figsize=(12, 5))
plt.scatter(data['Decimal_Year'], data['illum'], alpha=0.5)
plt.xlabel('Decimal Year')
plt.ylabel('Illumination Fraction')
plt.title('BBSO Illumination Fraction vs Decimal Year for 1998-2017')
plt.show()

# %%
