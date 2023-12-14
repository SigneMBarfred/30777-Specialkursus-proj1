
import pandas as pd
import numpy as np
import csv
import math

import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import stats

from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
import astropy.units as u

# Load IRAC data
IRAC_data = pd.read_csv('C:/Users/Nikolaj Lange Dons/OneDrive - Danmarks Tekniske Universitet/Dokumenter/5 semester/Space special kursus/IRAC_source_cat_3c239_400arcsec_radius.csv')

m1 = np.isnan(IRAC_data['i1_f_ap1']) & np.isnan(IRAC_data['i2_f_ap1'])
IRAC = IRAC_data[~m1]

# Constants Vega mag to flux
Fv0_W1 = 280.9
Fv0_W2 = 179.7
Fv0_W3 = 115
Fv0_W4 = 64.9

ch1 = 2.79
ch2 = 3.26
ch3 = 3.73
ch4 = 4.40

# Unit definition mJy
IRAC['i1_f_ap1'] = -1.085736205*np.log(IRAC['i1_f_ap1']/Fv0_W1) + ch1
IRAC['i2_f_ap1'] = -1.085736205*np.log(IRAC['i2_f_ap1']/Fv0_W2) + ch2
IRAC['i3_f_ap1'] = -1.085736205*np.log(IRAC['i3_f_ap1']/Fv0_W3) + ch3
IRAC['i4_f_ap1'] = -1.085736205*np.log(IRAC['i4_f_ap1']/Fv0_W4) + ch4

# Unit definition mJy - error
IRAC['i1_df_ap1'] = -1.085736205*np.log(IRAC['i1_df_ap1']/Fv0_W1) + ch1
IRAC['i2_df_ap1'] = -1.085736205*np.log(IRAC['i2_df_ap1']/Fv0_W2) + ch2
IRAC['i3_df_ap1'] = -1.085736205*np.log(IRAC['i3_df_ap1']/Fv0_W3) + ch3
IRAC['i4_df_ap1'] = -1.085736205*np.log(IRAC['i4_df_ap1']/Fv0_W4) + ch4


# Unit definition mJy
IRAC['i1_f_ap1'] = IRAC['i1_f_ap1'] * 10**(-6)
IRAC['i2_f_ap1'] = IRAC['i2_f_ap1'] * 10**(-6)
IRAC['i3_f_ap1'] = IRAC['i3_f_ap1'] * 10**(-6)
IRAC['i4_f_ap1'] = IRAC['i4_f_ap1'] * 10**(-6)

# Unit definition mJy - error
IRAC['i1_df_ap1'] = IRAC['i1_df_ap1'] * 10**(-6)
IRAC['i2_df_ap1'] = IRAC['i2_df_ap1'] * 10**(-6)
IRAC['i3_df_ap1'] = IRAC['i3_df_ap1'] * 10**(-6)
IRAC['i4_df_ap1'] = IRAC['i4_df_ap1'] * 10**(-6)

# IRAC['i1_mag'] = -1.085736205*np.log(IRAC['i1_f_ap1']/Fv0_W1) + ch1
# IRAC['i2_mag'] = -1.085736205*np.log(IRAC['i2_f_ap1']/Fv0_W2) + ch2
# IRAC['i3_mag'] = -1.085736205*np.log(IRAC['i3_f_ap1']/Fv0_W3) + ch3
# IRAC['i4_mag'] = -1.085736205*np.log(IRAC['i4_f_ap1']/Fv0_W4) + ch4

# # Unit definition mJy - error
# IRAC['i1_df_mag'] = -1.085736205*np.log(IRAC['i1_df_ap1']/Fv0_W1) + ch1
# IRAC['i2_df_mag'] = -1.085736205*np.log(IRAC['i2_df_ap1']/Fv0_W2) + ch2
# IRAC['i3_df_mag'] = -1.085736205*np.log(IRAC['i3_df_ap1']/Fv0_W3) + ch3
# IRAC['i4_df_mag'] = -1.085736205*np.log(IRAC['i4_df_ap1']/Fv0_W4) + ch4

# RAGERS coordinates txt file
positions_file = 'C:/Users/Nikolaj Lange Dons/OneDrive - Danmarks Tekniske Universitet/Dokumenter/5 semester/Space special kursus/RAGERStxt.txt'
txt_coord = pd.read_csv(positions_file, sep=" ")
obj_id = txt_coord['obj_id']

# Define coordinates
RAGER_ra = txt_coord['ra']
RAGER_dec = txt_coord['dec']
IRAC_coord = IRAC[['ra', 'dec']]

# Input: IRAC dataframe - RAGER ra and dec in deg - Other source ra and dec in deg - radius in arcsec
def ang_dist_IRAC(df, RAGER_ra, RAGER_dec, source_ra, source_dec, radius):
    df['angular_distance'] = np.arccos(
        np.sin(RAGER_dec * np.pi / 180) * np.sin(source_dec * np.pi / 180) +
        np.cos(RAGER_dec * np.pi / 180) * np.cos(source_dec * np.pi / 180) *
        np.cos(RAGER_ra * np.pi / 180 - source_ra * np.pi / 180)
    ) * 180 / np.pi
    df['mask'] = df['angular_distance'] <= radius
    return df['mask'], df['angular_distance']

# Create n dataframes: Will contain data on sources within 15 arcsec of RAGERS
RAGERS_IRAC = [pd.DataFrame() for _dummy in range(len(txt_coord['ra']))]

# Adds columns with names mask1 - mask24 including boolean if source is within target radius
# Adds columns with names ang_dist1 - ang_dist24 including boolean if source is within target radius
for i in range(0,len(txt_coord['ra'])):
    mask_id = ('mask' + str(i+1))
    ang_dist_id = ('ang_distance' + str(i+1))
    IRAC[str(mask_id)] = ang_dist_IRAC(IRAC, RAGER_ra[i], RAGER_dec[i], IRAC_coord['ra'], IRAC_coord['dec'], 15.0/3600.0)[0]
    IRAC[str(ang_dist_id)] = ang_dist_IRAC(IRAC, RAGER_ra[i], RAGER_dec[i], IRAC_coord['ra'], IRAC_coord['dec'], 15.0/3600.0)[1]
    
    # Loop to create list of dataframes with sources within target radius 15 arcsec of RAGERS
    RAGERS_IRAC[i] = pd.DataFrame(IRAC[IRAC.columns[0:56]][IRAC[str(mask_id)]])
    RAGERS_IRAC[i]['ang_dist'] = IRAC[str(ang_dist_id)][IRAC[str(mask_id)]]


# Count sources around each RAGER for IRAC
n_source_IRAC = np.array(np.zeros(len(txt_coord['ra'])))

for i in range(len(txt_coord['ra'])):
    if not RAGERS_IRAC[i].empty:
        n_source_IRAC[i] = len(RAGERS_IRAC[i])
        print(f'RAGER{i + 1}: obj_id {obj_id[i]} has {n_source_IRAC[i]} sources within 15 arcseconds')
    else:
        n_source_IRAC[i] = 0
        print(f'RAGER{i + 1}: obj_id {obj_id[i]} has {n_source_IRAC[i]} sources within 15 arcseconds')

# Get source name for sources within 15 arcseconds for IRAC
for i in range(len(txt_coord['ra'])):
    if not RAGERS_IRAC[i].empty:
        IRAC_names = RAGERS_IRAC[i]['objid'].values
        print(f'RAGER{i + 1}: obj_id {obj_id[i]} is near {IRAC_names} ')

# Loop to calculate flux for IRAC sources around RAGERS
flux_values_IRAC = []

for i in range(len(txt_coord['ra'])):
    flux_values_IRAC.append([obj_id[i]] + RAGERS_IRAC[i][['i1_f_ap1', 'i2_f_ap1', 'i3_f_ap1', 'i4_f_ap1']].values.flatten().tolist())

# Save the flux values to a text file
output_file_IRAC = 'flux_values_IRAC.txt'
with open(output_file_IRAC, 'w') as f:
    f.write('obj_id\tflux_I1[r\"$\mu Jy$"]\tflux_I2[r\"$\mu Jy$"]\tflux_I3[r\"$\mu Jy$"]\tflux_I4[r\"$\mu Jy$"]\n')
    for entry in flux_values_IRAC:
        f.write('\t'.join(map(str, entry)) + '\n')

print(f'Flux values saved to {output_file_IRAC}')

# 3. 7, 11, 14, 16

###########################
### color-color diagram ###

wl_ch1 = 3.6*u.micron
wl_ch2 = 4.5*u.micron
wl_ch3 = 5.8*u.micron
wl_ch4 = 8*u.micron

# Precise wavelength [microns]
# ch1: 3.5686
# ch2: 4.5067
# ch3: 5.7788
# ch4: 7.9958

### IRAC ###
 
# Equation - log Flux Density (muJy - Micro Jansky)
IRAC['S1-S2'] = np.log10(IRAC['i1_f_ap1']/IRAC['i2_f_ap1']) # Ch3-Ch1 - color
IRAC['S2'] = np.log10(IRAC['i2_f_ap1'])                     # Ch2 - flux
IRAC['S2-S3'] = np.log10(IRAC['i2_f_ap1']/IRAC['i3_f_ap1']) # Ch4-Ch2 - color
# IRAC['m1-m2'] = IRAC['w1mpro'] - IRAC['w2mpro']
# IRAC['m2-m3'] = IRAC['w2mpro'] - IRAC['w3mpro']

# Equation - log10 Flux Density (muJy - Micro Jansky) ERROR
IRAC['S1-S2_err'] = np.log10(IRAC['i1_df_ap1']/IRAC['i2_df_ap1']) # Ch3-Ch1 - color (ERROR)
IRAC['S2_err'] = np.log10(IRAC['i2_df_ap1'])                      # Ch2 - flux (ERROR)
IRAC['S2-S3_err'] = np.log10(IRAC['i2_df_ap1']/IRAC['i3_df_ap1']) # Ch4-Ch2 - color (ERROR)
# IRAC['m1-m2_err'] = IRAC['w1sigmpro'] - IRAC['w2sigmpro']
# IRAC['m2-m3_err'] = IRAC['w2sigmpro'] - IRAC['w3sigmpro']

### RAGERS_IRAC ###
for i in range(0,len(txt_coord['ra'])):    
    # Equation 1 - log10 Flux Density (muJy - Micro Jansky) and color
    RAGERS_IRAC[i]['S1-S2'] = np.log10(RAGERS_IRAC[i]['i1_f_ap1']/RAGERS_IRAC[i]['i2_f_ap1'])
    RAGERS_IRAC[i]['S2'] = np.log10(RAGERS_IRAC[i]['i2_f_ap1'])
    RAGERS_IRAC[i]['S2-S3'] = np.log10(RAGERS_IRAC[i]['i2_f_ap1']/RAGERS_IRAC[i]['i3_f_ap1'])
    
    # Equation 1 - log10 Flux Density (muJy - Micro Jansky) and color (ERROR)
    RAGERS_IRAC[i]['S1-S2_err'] = np.log10(RAGERS_IRAC[i]['i1_df_ap1']/RAGERS_IRAC[i]['i2_df_ap1'])
    RAGERS_IRAC[i]['S2_err'] = np.log10(RAGERS_IRAC[i]['i2_df_ap1'])/np.log10(10)
    RAGERS_IRAC[i]['S2-S3_err'] = np.log10(RAGERS_IRAC[i]['i2_df_ap1']/RAGERS_IRAC[i]['i3_df_ap1'])

# Color arrays    
color_RAGERS_IRAC_x = np.array([RAGERS_IRAC[2]['S1-S2'].values[0],RAGERS_IRAC[2]['S1-S2'].values[1],
    RAGERS_IRAC[13]['S1-S2'].values[0], RAGERS_IRAC[15]['S1-S2'].values[0]])

color_RAGERS_IRAC_y = np.array([RAGERS_IRAC[2]['S2-S3'].values[0],RAGERS_IRAC[2]['S2-S3'].values[1],
    RAGERS_IRAC[13]['S2-S3'].values[0], RAGERS_IRAC[15]['S2-S3'].values[0]])

flux_RAGERS_IRAC_y = np.array([RAGERS_IRAC[2]['S2'].values[0],RAGERS_IRAC[2]['S2'].values[1],
    RAGERS_IRAC[13]['S2'].values[0], RAGERS_IRAC[15]['S2'].values[0]])

# 3. 7, 11, 14, 16
# Color arrays - error
color_RAGERS_IRAC_err_x = np.array([RAGERS_IRAC[2]['S1-S2_err'].values[0],RAGERS_IRAC[2]['S1-S2_err'].values[1],
    RAGERS_IRAC[13]['S1-S2_err'].values[0], RAGERS_IRAC[15]['S1-S2_err'].values[0]])

color_RAGERS_IRAC_err_y = np.array([RAGERS_IRAC[2]['S2-S3_err'].values[0],RAGERS_IRAC[2]['S2-S3_err'].values[1],
    RAGERS_IRAC[13]['S2-S3_err'].values[0], RAGERS_IRAC[15]['S2-S3_err'].values[0]])

flux_RAGERS_IRAC_err_y = np.array([RAGERS_IRAC[2]['S2_err'].values[0],RAGERS_IRAC[2]['S2_err'].values[1],
    RAGERS_IRAC[13]['S2_err'].values[0], RAGERS_IRAC[15]['S2_err'].values[0]])


# Arrays with figure info about color, marker type and label
color = ['black', 'darkred', 'darkorange', 'steelblue']
marker = ['o', 'o', 'o', 'o']
label = ['3(1)','3(2)','14','16']


# x: Ch3-Ch1 color, y: Ch4-Ch2 color
fig, ax = plt.subplots()
hex1 = ax.hexbin(IRAC['S1-S2'], IRAC['S2'], vmax = 1, cmap = "binary", mincnt = 0, gridsize=(173,100))

for i in [0,1,2,3]:
    plt.scatter(color_RAGERS_IRAC_x[i], flux_RAGERS_IRAC_y[i], marker=marker[i], color=color[i], label=label[i])
    # plt.errorbar(color_RAGERS_IRAC_x[i], color_RAGERS_IRAC_y[i], yerr=color_RAGERS_IRAC_err_y[i], fmt=marker[i], color=color[i])
    # plt.errorbar(color_RAGERS_IRAC_x[i], color_RAGERS_IRAC_y[i], xerr=color_RAGERS_IRAC_err_x[i], fmt=marker[i], color=color[i])
    
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# plt.plot([0.22, 0.22], [-5.3, -3.7], color= 'r' , linestyle='dashed')   

plt.title('IRAC color diagram', fontsize=13)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="RAGER obj_id") 
plt.xlabel("Ch1 - Ch2", fontsize=10)
plt.ylabel("Ch2", fontsize=10)
# plt.xlim(-1, 3)
# plt.ylim(0,6)

plt.show()

# color cut

### color - color cuts ###


color_mask = IRAC['S1-S2'] > 0.22
n_total = np.sum(color_mask)
color_ratio = n_total/len(IRAC['S1-S2'])

SMG_colors = [['SMG_color'+ str(_dummy)] for _dummy in range(1, 25)]

for i in range(0,len(txt_coord['ra'])):
    mask_id = ('mask' + str(i+1))
    mask_id1 = ('SMG_mask' + str(i+1))
    IRAC[str(mask_id1)] = color_mask & IRAC[str(mask_id)]
    
    SMG_colors[i] = pd.DataFrame(IRAC[IRAC.columns[0:32]][IRAC[str(mask_id1)]])


for i in range(len(txt_coord['ra'])):
    if any(SMG_colors[i]['ra']):
        print([i, SMG_colors[i]['designation']])



### Histogram ###

# Hist1 for Ch3-Ch1
plt.figure(figsize=(12,8), edgecolor='blue')
plt.hist(IRAC['S1-S2'], bins = 20)
# plt.ylim(0,30)
plt.title('Hist for IRAC-color Ch1-Ch2', fontsize = 16)
plt.xlabel('Ch1 - Ch2 flux-color', fontsize = 14)
plt.ylabel('Number of occourence', fontsize = 14)
plt.show()

# Hist2 for Ch4-Ch2
plt.figure(figsize=(12, 8), edgecolor='blue')
plt.hist(IRAC['S2'], bins = 20 )
# plt.ylim(0,40)
plt.title('Hist for IRAC-flux Ch2', fontsize = 16)
plt.xlabel('IRAC-flux Ch2', fontsize = 14)
plt.ylabel('Number of occourence', fontsize = 14)
plt.show()


### RAGERS calling commands ###

# RAGERS[RAGER number with base zero][string column name]
# To call multiple columns:
    # RAGERS[15][['ang_dist', 'flux_W1']]

# Get only value from dataframe cell:
# RAGERS[i]['designation'].values[0]

############################
### Unused lines of code ###

# # Define coordinates as skycoord
# c1 = SkyCoord(ra=x_coord*u.degree, dec=y_coord*u.degree, distance=1500.3*u.pc) # RAGERS
# c2 = SkyCoord(ra=WISE_coord['ra']*u.degree, dec=WISE_coord['dec']*u.degree, distance=WISE['dist']*u.pc) # Other source

# Function to ceate new columns with distance to RAGERS
# def dist_to_RAGER(RAGERS_ra, RAGERS_dec, RAGERS_dist, WISE_ra, WISE_dec, WISE_dist):
    
#     c1 = SkyCoord(ra=RAGERS_ra*u.degree, dec=RAGERS_dec*u.degree, distance=RAGERS_dist*u.pc) # RAGERS
#     c2 = SkyCoord(ra=WISE_ra*u.degree, dec=WISE_dec*u.degree, distance=WISE_dist*u.pc)

#     # Obtain astropy's distance between c1 & c2 coords.
#     radial_dist = c1.separation_3d(c2)
#     return radial_dist

# Calculate distance to RAGER sources with Hubble equation
# def Hubble_dist(z):
#     c = 299792458*u.m*u.s**(-1)
#     v = c * z
    
#     H_0 = 70500*u.m*u.Mpc**(-1)*u.s**(-1)
#     dist = v/H_0

#     return dist

# x: Ch3-Ch1 color, y: Ch4-Ch2 color
fig, ax = plt.subplots()
plt.scatter(color_RAGERS_IRAC_x[0], color_RAGERS_IRAC_y[0], marker=marker[0], color=color[0], label=label[0])

plt.errorbar(color_RAGERS_IRAC_x[0], color_RAGERS_IRAC_y[0], yerr=color_RAGERS_IRAC_err_y[0])
# plt.errorbar(color_RAGERS_IRAC_x[0], color_RAGERS_IRAC_y[0], xerr=color_RAGERS_IRAC_err_x[0])

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set(title='Galaxy color-color diagram')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="RAGER obj_id") 
plt.xlabel("log $S_{5.8}$/$S_{3.6}$")
plt.ylabel("log $S_{8.0}$/$S_{4.5}$")
plt.show()

# loop color error
for i in range(13):
    print([i+1, color_RAGERS_IRAC_x[i], color_RAGERS_IRAC_y[i]])

a = color_RAGERS_IRAC_x
b = color_RAGERS_IRAC_y
plt.scatter(a, b)

c1 = color_RAGERS_IRAC_err_x
c2 = color_RAGERS_IRAC_err_y[0]
 
plt.errorbar(a, b, xerr=c1, fmt="o")
plt.errorbar(a, b, yerr=c2, fmt="o")

plt.show()




plt.figure()

a = color_RAGERS_IRAC_x
b = color_RAGERS_IRAC_y
plt.scatter(a, b)
 
c = color_RAGERS_IRAC_err_x
d = color_RAGERS_IRAC_err_y

plt.errorbar(a, b, xerr=c, yerr=d, fmt="o", color="r")
 
plt.show()


[2, 6, 10, 13, 15]

for i in [11]:
    print('W1')
    print(RAGERS_IRAC[i]['i1_f_ap1'])
    print('W2')
    print(RAGERS_IRAC[i]['i1_f_ap1'])
    print('W3')
    print(RAGERS_IRAC[i]['i1_f_ap1'])
    print('W4')
    print(RAGERS_IRAC[i]['i1_f_ap1'])





