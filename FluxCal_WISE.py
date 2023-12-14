# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 21:56:42 2023

@author: Nikolaj Lange Dons
"""

################
### Get Flux ###
################

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

# load Scuba data
# hdu_scuba = fits.open('C:/Users/Nikolaj Lange Dons/OneDrive - Danmarks Tekniske Universitet/Dokumenter/5 semester/Space special kursus/Saved files/cutout1.0.fits')
hdu_scuba = fits.open('C:/Users/Nikolaj Lange Dons/OneDrive - Danmarks Tekniske Universitet/Dokumenter/5 semester/Space special kursus/3C239_850_mf_cal_crop_snr.fits')
scuba850_data_cube = hdu_scuba[0].data  # This is a 3D cube
hdr_scuba = hdu_scuba[0].header
wcs_scuba = WCS(hdr_scuba)
df_scuba = pd.DataFrame(scuba850_data_cube[0,:,:]) # Pandas dataframe

# RAGERS koordinater txt fil
positions_file = 'C:/Users/Nikolaj Lange Dons/OneDrive - Danmarks Tekniske Universitet/Dokumenter/5 semester/Space special kursus/RAGERStxt.txt'
txt_coord = pd.read_csv(positions_file, sep=" ") # Pandas dataframe
obj_id = txt_coord['obj_id']

# WISE data
WISE = pd.read_csv('C:/Users/Nikolaj Lange Dons/OneDrive - Danmarks Tekniske Universitet/Dokumenter/5 semester/Space special kursus/Saved files/table_irsa_catalog_search_results.csv')
WISE_header = WISE.columns

##########################################
### calculate flux and AB mag for WISE ###

# Constants Vega mag to flux
Fv0_W1 = 309.540
Fv0_W2 = 171.757
Fv0_W3 = 31.678
Fv0_W4 = 8.363

### WISE ###
# Equation 1 - Vega Magnitudes to Flux Density muJy - Micro Jansky
WISE['flux_W1'] = Fv0_W1 * 10**(-WISE['w1mpro']/2.5)*10**6
WISE['flux_W2'] = Fv0_W2 * 10**(-WISE['w2mpro']/2.5)*10**6
WISE['flux_W3'] = Fv0_W3 * 10**(-WISE['w3mpro']/2.5)*10**6
WISE['flux_W4'] = Fv0_W4 * 10**(-WISE['w4mpro']/2.5)*10**6

# Equation 1 - Vega Magnitudes uncertainty to Flux Density uncertainty (muJy - Micro Jansky)
WISE['flux_error_W1'] = Fv0_W1 * 10**(-WISE['w1sigmpro']/2.5)*10**6
WISE['flux_error_W2'] = Fv0_W2 * 10**(-WISE['w2sigmpro']/2.5)*10**6
WISE['flux_error_W3'] = Fv0_W3 * 10**(-WISE['w3sigmpro']/2.5)*10**6
WISE['flux_error_W4'] = Fv0_W4 * 10**(-WISE['w4sigmpro']/2.5)*10**6


#####################
# Find target sources around RAGERS

# Input: WISE dataframe - RAGER ra and dec in deg - Other source ra and dec in deg - radius in arcsec
def ang_dist(df, RAGER_ra, RAGER_dec, source_ra, source_dec, radius):
    df['angular_distance'] = np.arccos(np.sin(RAGER_dec*np.pi/180)*np.sin(source_dec*np.pi/180)+np.cos(RAGER_dec*np.pi/180)*np.cos(source_dec*np.pi/180)*np.cos(RAGER_ra*np.pi/180 - source_ra*np.pi/180))*180/np.pi
    df['mask'] = df['angular_distance'] <= radius
    return df['mask'], df['angular_distance']

# Define coordinates
RAGER_ra = txt_coord['ra'] # RAGERS x_coord
RAGER_dec = txt_coord['dec'] # RAGERS y_coord
WISE_coord = WISE[['ra','dec']]

# Create n dataframes: Will contain data on sources within 15 arcsec of RAGERS
RAGERS = [['RAGER'+ str(_dummy)] for _dummy in range(1, 25)]

# Adds columns with names mask1 - mask24 including boolean if source is within target radius
# Adds columns with names ang_dist1 - ang_dist24 including boolean if source is within target radius
for i in range(0,len(txt_coord['ra'])):
    mask_id = ('mask' + str(i+1))
    ang_dist_id = ('ang_dist' + str(i+1))
    WISE[str(mask_id)] = ang_dist(WISE, RAGER_ra[i], RAGER_dec[i], WISE_coord['ra'], WISE_coord['dec'], 15.0/3600.0)[0]
    WISE[str(ang_dist_id)] = ang_dist(WISE, RAGER_ra[i], RAGER_dec[i], WISE_coord['ra'], WISE_coord['dec'], 15.0/3600.0)[1]
    
    # Loop to create list of dataframes with sources within target radius 15 arcsec of RAGERS
    RAGERS[i] = pd.DataFrame(WISE[WISE.columns[0:32]][WISE[str(mask_id)]])
    RAGERS[i]['ang_dist'] = WISE[str(ang_dist_id)][WISE[str(mask_id)]]

# Count sources around each RAGER
n_source = np.array(np.zeros(len(txt_coord['ra'])))

for i in range(0,len(txt_coord['ra'])):
    
    if any(RAGERS[i]['w1mpro']):
        n_source[i] = np.array([len(RAGERS[i]['w1mpro'])])
        print(f'RAGER{i+1}: obj_id {obj_id[i]} has {n_source[i]} sources within 15 arcseconds')
    
    else:
        n_source[i] = 0
        print(f'RAGER{i+1}: obj_id {obj_id[i]} has {n_source[i]} sources within 15 arcseconds')

# Get source name for sources within 15 arcseconds
for i in range(0,len(txt_coord['ra'])):
    
    if any(RAGERS[i]['w1mpro']):
        WISE_names = RAGERS[i]['designation'].values[:]
        print(f'RAGER{i+1}: obj_id {obj_id[i]} is near {WISE_names} ')

# SCUBA/RAGER
# z = 1.781


### RAGERS ###
for i in range(0,len(txt_coord['ra'])):
    # Define vega magnitude for all chanels
    m_vega_W1 = RAGERS[i]['w1mpro']
    m_vega_W2 = RAGERS[i]['w2mpro']
    m_vega_W3 = RAGERS[i]['w3mpro']
    m_vega_W4 = RAGERS[i]['w4mpro']

    # Flux uncerntanty in Vega Magnitude  units for all chanels
    m_vega_error_W1 = RAGERS[i]['w1sigmpro']
    m_vega_error_W2 = RAGERS[i]['w2sigmpro']
    m_vega_error_W3 = RAGERS[i]['w3sigmpro']
    m_vega_error_W4 = RAGERS[i]['w4sigmpro']

    # Equation 1 - Vega Magnitudes to Flux Density
    RAGERS[i]['flux_W1'] = Fv0_W1 * 10**(-m_vega_W1/2.5)*10**6
    RAGERS[i]['flux_W2'] = Fv0_W2 * 10**(-m_vega_W2/2.5)*10**6
    RAGERS[i]['flux_W3'] = Fv0_W3 * 10**(-m_vega_W3/2.5)*10**6
    RAGERS[i]['flux_W4'] = Fv0_W4 * 10**(-m_vega_W4/2.5)*10**6

    # Equation 1 - Vega Magnitudes uncertainty to Flux Density uncertainty
    RAGERS[i]['flux_error_W1'] = Fv0_W1 * 10**(-m_vega_error_W1/2.5)*10**6
    RAGERS[i]['flux_error_W2'] = Fv0_W2 * 10**(-m_vega_error_W2/2.5)*10**6
    RAGERS[i]['flux_error_W3'] = Fv0_W3 * 10**(-m_vega_error_W3/2.5)*10**6
    RAGERS[i]['flux_error_W4'] = Fv0_W4 * 10**(-m_vega_error_W4/2.5)*10**6

for i in range(23):
    if any(RAGERS[i]['flux_error_W1']):
        print([i, RAGERS[i]['flux_error_W1'].values])

for i in range(23):
    if any(RAGERS[i]['flux_error_W2']):
        print([i, RAGERS[i]['flux_error_W2'].values])

for i in range(23):
    if any(RAGERS[i]['flux_error_W3']):
        print([i, RAGERS[i]['flux_error_W3'].values])

for i in range(23):
    if any(RAGERS[i]['flux_error_W4']):
        print([i, RAGERS[i]['flux_error_W4'].values])
    

############################
### Galaxy color-diagram ###

# wavelength [microns]
wl_ch1 = 3.4*u.micron
wl_ch2 = 4.6*u.micron
wl_ch3 = 11.5*u.micron
wl_ch4 = 22*u.micron

### WISE ###
 
# Equation - log Flux Density (muJy - Micro Jansky)
WISE['S1-S2'] = np.log10(WISE['flux_W1']/WISE['flux_W2']) # Ch2/Ch1 - color
WISE['S3'] = np.log10(WISE['flux_W3'])                    # Ch3 - flux
WISE['S2-S3'] = np.log10(WISE['flux_W2']/WISE['flux_W3']) # Ch3/Ch2 - color
WISE['m1-m2'] = WISE['w1mpro'] - WISE['w2mpro']
WISE['m2-m3'] = WISE['w2mpro'] - WISE['w3mpro']

# Equation - log Flux Density (muJy - Micro Jansky) ERROR
WISE['S1-S2_err'] = np.log10(WISE['flux_error_W1']/WISE['flux_error_W2']) # Ch2/Ch1 - color (ERROR)
WISE['S3_err'] = np.log10(WISE['flux_error_W3'])                          # Ch3 - flux (ERROR)
WISE['S2-S3_err'] = np.log10(WISE['flux_error_W2']/WISE['flux_error_W3']) # Ch3/Ch2 - color (ERROR)
WISE['m1-m2_err'] = WISE['w1sigmpro'] - WISE['w2sigmpro']
WISE['m2-m3_err'] = WISE['w2sigmpro'] - WISE['w3sigmpro']

### RAGERS ###
for i in range(0,len(txt_coord['ra'])):    
    # Equation 1 - log10 Flux Density (muJy - Micro Jansky) and color
    RAGERS[i]['S1-S2'] = np.log10(RAGERS[i]['flux_W1']/RAGERS[i]['flux_W2'])
    RAGERS[i]['S3'] = np.log10(RAGERS[i]['flux_W3'])
    RAGERS[i]['S2-S3'] = np.log10(RAGERS[i]['flux_W2']/RAGERS[i]['flux_W3'])
    RAGERS[i]['m1-m2'] = RAGERS[i]['w1mpro'] - RAGERS[i]['w2mpro']
    RAGERS[i]['m2-m3'] = RAGERS[i]['w2mpro'] - RAGERS[i]['w3mpro']
    
    # Equation 1 - log10 Flux Density (muJy - Micro Jansky) and color (ERROR)
    RAGERS[i]['S1-S2_err'] = np.log10(RAGERS[i]['flux_error_W1']/RAGERS[i]['flux_error_W2'])
    RAGERS[i]['S3_err'] = np.log10(RAGERS[i]['flux_error_W3'])
    RAGERS[i]['S2-S3_err'] = np.log10(RAGERS[i]['flux_error_W2']/RAGERS[i]['flux_error_W3'])
    RAGERS[i]['m1-m2_err'] = RAGERS[i]['w1sigmpro'] - RAGERS[i]['w2sigmpro']
    RAGERS[i]['m2-m3_err'] = RAGERS[i]['w2sigmpro'] - RAGERS[i]['w3sigmpro']

# Color arrays    
color_RAGERS_y = np.array([RAGERS[0]['S1-S2'].values[0],RAGERS[2]['S1-S2'].values[0],RAGERS[5]['S1-S2'].values[0],RAGERS[6]['S1-S2'].values[0],
    RAGERS[8]['S1-S2'].values[0], RAGERS[11]['S1-S2'].values[0],RAGERS[12]['S1-S2'].values[0],RAGERS[13]['S1-S2'].values[0],
    RAGERS[14]['S1-S2'].values[0], RAGERS[15]['S1-S2'].values[0], RAGERS[15]['S1-S2'].values[1], RAGERS[17]['S1-S2'].values[0], RAGERS[21]['S1-S2'].values[0]])

color_RAGERS_x = np.array([RAGERS[0]['S2-S3'].values[0],RAGERS[2]['S2-S3'].values[0],RAGERS[5]['S2-S3'].values[0],RAGERS[6]['S2-S3'].values[0],
    RAGERS[8]['S2-S3'].values[0], RAGERS[11]['S2-S3'].values[0],RAGERS[12]['S2-S3'].values[0],RAGERS[13]['S2-S3'].values[0],
    RAGERS[14]['S2-S3'].values[0], RAGERS[15]['S2-S3'].values[0], RAGERS[15]['S2-S3'].values[1], RAGERS[17]['S2-S3'].values[0], RAGERS[21]['S2-S3'].values[0]])

flux_RAGERS_y = np.array([RAGERS[0]['S3'].values[0],RAGERS[2]['S3'].values[0],RAGERS[5]['S3'].values[0],RAGERS[6]['S3'].values[0],
    RAGERS[8]['S3'].values[0], RAGERS[11]['S3'].values[0],RAGERS[12]['S3'].values[0],RAGERS[13]['S3'].values[0],
    RAGERS[14]['S3'].values[0], RAGERS[15]['S3'].values[0], RAGERS[15]['S3'].values[1], RAGERS[17]['S3'].values[0], RAGERS[21]['S3'].values[0]])

color_mag_x = np.array([RAGERS[0]['m2-m3'].values[0],RAGERS[2]['m2-m3'].values[0],RAGERS[5]['m2-m3'].values[0],RAGERS[6]['m2-m3'].values[0],
    RAGERS[8]['m2-m3'].values[0], RAGERS[11]['m2-m3'].values[0],RAGERS[12]['m2-m3'].values[0],RAGERS[13]['m2-m3'].values[0],
    RAGERS[14]['m2-m3'].values[0], RAGERS[15]['m2-m3'].values[0], RAGERS[15]['m2-m3'].values[1], RAGERS[17]['m2-m3'].values[0], RAGERS[21]['m2-m3'].values[0]])

color_mag_y = np.array([RAGERS[0]['m1-m2'].values[0],RAGERS[2]['m1-m2'].values[0],RAGERS[5]['m1-m2'].values[0],RAGERS[6]['m1-m2'].values[0],
    RAGERS[8]['m1-m2'].values[0], RAGERS[11]['m1-m2'].values[0],RAGERS[12]['m1-m2'].values[0],RAGERS[13]['m1-m2'].values[0],
    RAGERS[14]['m1-m2'].values[0], RAGERS[15]['m1-m2'].values[0], RAGERS[15]['m1-m2'].values[1], RAGERS[17]['m1-m2'].values[0], RAGERS[21]['m1-m2'].values[0]])

# Color arrays - error
color_RAGERS_err_y = np.array([RAGERS[0]['S1-S2_err'].values[0],RAGERS[2]['S1-S2_err'].values[0],RAGERS[5]['S1-S2_err'].values[0],RAGERS[6]['S1-S2_err'].values[0],
    RAGERS[8]['S1-S2_err'].values[0], RAGERS[11]['S1-S2_err'].values[0],RAGERS[12]['S1-S2_err'].values[0],RAGERS[13]['S1-S2_err'].values[0],
    RAGERS[14]['S1-S2_err'].values[0], RAGERS[15]['S1-S2_err'].values[0], RAGERS[15]['S1-S2_err'].values[1], RAGERS[17]['S1-S2_err'].values[0], RAGERS[21]['S1-S2_err'].values[0]])

color_RAGERS_err_x = np.array([RAGERS[0]['S2-S3_err'].values[0],RAGERS[2]['S2-S3_err'].values[0],RAGERS[5]['S2-S3_err'].values[0],RAGERS[6]['S2-S3_err'].values[0],
    RAGERS[8]['S2-S3_err'].values[0], RAGERS[11]['S2-S3_err'].values[0],RAGERS[12]['S2-S3_err'].values[0],RAGERS[13]['S2-S3_err'].values[0],
    RAGERS[14]['S2-S3_err'].values[0], RAGERS[15]['S2-S3_err'].values[0], RAGERS[15]['S2-S3_err'].values[1], RAGERS[17]['S2-S3_err'].values[0], RAGERS[21]['S2-S3_err'].values[0]])

flux_RAGERS_err_y = np.array([RAGERS[0]['S3_err'].values[0],RAGERS[2]['S3_err'].values[0],RAGERS[5]['S3_err'].values[0],RAGERS[6]['S3_err'].values[0],
    RAGERS[8]['S3_err'].values[0], RAGERS[11]['S3_err'].values[0],RAGERS[12]['S3_err'].values[0],RAGERS[13]['S3_err'].values[0],
    RAGERS[14]['S3_err'].values[0], RAGERS[15]['S3_err'].values[0], RAGERS[15]['S3_err'].values[1], RAGERS[17]['S3_err'].values[0], RAGERS[21]['S3_err'].values[0]])

# Arrays with figure info about color, marker type and label
colors = ['black', 'grey','lightcoral', 'darkred', 'darkorange', 'khaki', 'forestgreen', 'lightseagreen', 'steelblue', 'navy', 'red', 'pink', 'coral']
markers = ['o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o']
labels = ['1','3','6','7','9','12','13','14','15','16(1)','16(2)','20','30']


# Combine the two previous plots
df = pd.DataFrame({'W2-W3':color_mag_x, 'W1-W2':color_mag_y,'label': ('1','3','6','7','9','12','13','14','15','16(1)','16(2)','20','30'),'color': sns.color_palette("bright", 13)})

fig, ax = plt.subplots()
# hex1 = ax.hexbin(WISE['m2-m3'], WISE['m1-m2'], vmax = 1, cmap = "binary", mincnt = 0, gridsize=(173,100), label='WISE')
scat1 = sns.scatterplot(x='W2-W3',y='W1-W2',hue = 'label', palette= sns.color_palette("bright", 13), data = df)
ax.set(title='Vega magntude color-color diagram [WISE]')
legend1 = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="RAGER obj_id")
# plt.xlim(2, 5)
# plt.ylim(-0.5, 1.6)
plt.show()

# x: Ch3-Ch1 color, y: Ch3 flux
fig, ax = plt.subplots()
hex1 = ax.hexbin(WISE['m2-m3'], WISE['m1-m2'], vmax = 1, cmap = "binary", mincnt = 0, gridsize=(173,100))

for i in range(13):
    plt.scatter(color_RAGERS_x[i], flux_RAGERS_y[i], marker=markers[i], color=colors[i], label=labels[i])
    
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set(title='color-flux diagram')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="RAGER obj_id") 
plt.xlabel("W2 - W3", fontsize=11)
plt.ylabel("W1 - W2", fontsize=11)
plt.show()

# x: Ch2-Ch3 color, y: Ch1-Ch2 color
fig, ax = plt.subplots()
hex1 = ax.hexbin(WISE['m2-m3'], WISE['m1-m2'], vmax = 1, cmap = "binary", mincnt = 0, gridsize=(173,100))

plt.scatter(color_mag_x[0], color_mag_y[0], marker=markers[0], color=colors[0], label=labels[0])
plt.scatter(color_mag_x[1], color_mag_y[1], marker=markers[1], color=colors[1], label=labels[1])
plt.scatter(color_mag_x[2], color_mag_y[2], marker=markers[2], color=colors[2], label=labels[2])
plt.scatter(color_mag_x[3], color_mag_y[3], marker=markers[3], color=colors[3], label=labels[3])
plt.scatter(color_mag_x[4], color_mag_y[4], marker=markers[4], color=colors[4], label=labels[4])
plt.scatter(color_mag_x[5], color_mag_y[5], marker=markers[5], color=colors[5], label=labels[5])
plt.scatter(color_mag_x[6], color_mag_y[6], marker=markers[6], color=colors[6], label=labels[6])
plt.scatter(color_mag_x[7], color_mag_y[7], marker=markers[7], color=colors[7], label=labels[7])
plt.scatter(color_mag_x[8], color_mag_y[8], marker=markers[8], color=colors[8], label=labels[8])
plt.scatter(color_mag_x[9], color_mag_y[9], marker=markers[9], color=colors[9], label=labels[9])
plt.scatter(color_mag_x[10], color_mag_y[10], marker=markers[10], color=colors[10], label=labels[10])
plt.scatter(color_mag_x[11], color_mag_y[11], marker=markers[11], color=colors[11], label=labels[11])
plt.scatter(color_mag_x[12], color_mag_y[12], marker=markers[12], color=colors[12], label=labels[12])

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('Vega magnitude color-color diagram [WISE]', fontsize=13)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="RAGER obj_id") 
plt.xlabel("W2 - W3", fontsize=11)
plt.ylabel("W1 - W2", fontsize=11)
plt.show()


### VEGA MAGNITUDE ###
df1 = pd.DataFrame({'x_axis':color_mag_x, 'y_axis':color_mag_y,'label': ('1','3','6','7','9','12','13','14','15','16(1)','16(2)','20','30'),'color': sns.color_palette("bright", 13)})

# x: Ch2-Ch3 color, y: Ch1-Ch2 color
fig, ax = plt.subplots()
hex1 = ax.hexbin(WISE['m2-m3'], WISE['m1-m2'], vmax = 1, cmap = "binary", mincnt = 0, gridsize=(173,100))

for i in range(13):
    plt.scatter(color_mag_x[i], color_mag_y[i], marker=markers[i], color=colors[i], label=labels[i])
    # plt.errorbar(color_RAGERS_x[i], color_RAGERS_y[i], yerr=color_RAGERS_err_y[i], fmt=marker[i], color=color[i])
    # plt.errorbar(color_RAGERS_x[i], color_RAGERS_y[i], xerr=color_RAGERS_err_x[i], fmt=marker[i], color=color[i])
    
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('color-color diagram [WISE]', fontsize=13)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="RAGER obj_id") 
plt.xlabel("W2 - W3", fontsize=11)
plt.ylabel("W1 - W2", fontsize=11)
plt.show()

fig, ax = plt.subplots()
# hex1 = ax.hexbin(WISE['m2-m3'], WISE['m1-m2'], vmax = 1, cmap = "binary", mincnt = 0, gridsize=(173,100))

for i in range(13):
    plt.scatter(color_mag_x[i], color_mag_y[i], marker=markers[i], color=colors[i], label=labels[i])
    print([i, colors[i], labels[i]])
    
    
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('color-color diagram [WISE]', fontsize=13)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="RAGER obj_id") 
plt.xlabel("W2 - W3", fontsize=11)
plt.ylabel("W1 - W2", fontsize=11)
plt.xlim(0, 5)
plt.ylim(-1, 2)
plt.show()


### color - color cuts ###

color_mask1 = WISE['m1-m2'] > WISE['m2-m3'] + 1.25 * WISE['m1-m2']
color_mask2 = WISE['m2-m3'] > 5.3
color_mask = color_mask1 & color_mask2

n_color = np.sum(color_mask)
color_ratio = n_color/len(WISE['m1-m2'])

SMG_colors = [['SMG_color'+ str(_dummy)] for _dummy in range(1, 25)]

for i in range(0,len(txt_coord['ra'])):
    mask_id = ('mask' + str(i+1))
    mask_id1 = ('SMG_mask' + str(i+1))
    WISE[str(mask_id1)] = color_mask & WISE[str(mask_id)]
    
    SMG_colors[i] = pd.DataFrame(WISE[WISE.columns[0:32]][WISE[str(mask_id1)]])    


for i in range(len(txt_coord['ra'])):
    if any(SMG_colors[i]['ra']):
        print([i, SMG_colors[i]['designation']])

### color-color diagram with color cuts ###

def graph(formula, x_range):  
    x = np.array(x_range)  
    y = formula(x)  # <- note now we're calling the function 'formula' with x
    plt.plot(x, y, color= 'r' , linestyle='dashed')   

def my_formula(x):
    return (7-x)/1.25

fig, ax = plt.subplots()
hex1 = ax.hexbin(WISE['m2-m3'], WISE['m1-m2'], vmax = 1, cmap = "binary", mincnt = 0, gridsize=(173,100))

for i in range(13):
    plt.scatter(color_mag_x[i], color_mag_y[i], marker=markers[i], color=colors[i], label=labels[i])
    # plt.errorbar(color_RAGERS_x[i], color_RAGERS_y[i], yerr=color_RAGERS_err_y[i], fmt=marker[i], color=color[i])
    # plt.errorbar(color_RAGERS_x[i], color_RAGERS_y[i], xerr=color_RAGERS_err_x[i], fmt=marker[i], color=color[i])

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))    
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

graph(my_formula, range(4, 7))

plt.title('color-color diagram [WISE]', fontsize=13)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="RAGER obj_id") 
plt.xlabel("W2 - W3", fontsize=11)
plt.ylabel("W1 - W2", fontsize=11)
plt.show()



### Histogram ###

# Hist1 for Ch1-Ch2
plt.figure(figsize=(12,8), edgecolor='blue')
plt.hist(WISE['m1-m2'], bins = 100)
# plt.ylim(0,30)
plt.title('Hist for WISE-color W1-W1', fontsize = 16)
plt.xlabel('W1-W2 color-color', fontsize = 14)
plt.ylabel('Number of occourence', fontsize = 14)
plt.show()

# Hist2 for Ch2-Ch3
plt.figure(figsize=(12,8), edgecolor='blue')
plt.hist(WISE['m2-m3'], bins = 100)
# plt.ylim(0,40)
plt.title('Hist for WISE-color W2-W3', fontsize = 16)
plt.xlabel('WISE-color W2-W3', fontsize = 14)
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


# # color-color diagram of WISE sources within 15 arcsec of a RAGERS
# plt.figure()
# for i in range(13):
#     plt.scatter(color_RAGERS_x[i], flux_RAGERS_y[i], marker=marker[i], color=color[i], label=label[i])
    
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.show() 

# # color-color diagram of WISE sources within 15 arcsec of a RAGERS
# plt.figure()
# navn = plt.scatter(color_RAGERS_x, flux_RAGERS_y, c=range(13), cmap='summer')
# plt.scatter(color_RAGERS_x[6], flux_RAGERS_y[6], marker='x', color='black')
# plt.title("Galaxy color-color diagram")
# plt.xlabel("log $S_{5.8}$/$S_{3.6}$") #defines x-axis label
# plt.ylabel("log $S_{5.8}$ ($\mu Jy$)") #defines y-axis label
# # plt.xlim(0, 2.5)
# # plt.ylim(-3.7, -3.2)
# plt.legend(handles=navn.legend_elements()[0], labels=('1','3','6','7','9','12','13','14','15','16(1)','16(2)','20','30') , loc='center left', bbox_to_anchor=(1, 0.5))
# plt.show

# loop color error


# # Hexbin af WISE - all sky catalouge - 600 arcsec radius af 3C239
# plt.figure()
# plt.hexbin(WISE['S3-S1'], WISE['S3'], vmax = 1, cmap = "binary", mincnt = 0, gridsize=(173,100))
# # set_xlim(0,2.5)
# # set_ylim(0,2.5)
# plt.ylabel("$S_{3}$", fontsize=12)
# plt.xlabel("$S_{3-1}$")
# plt.title('Color-Color diagram (WISE)')
# plt.show()


# for i in range(13):
#     print([i+1, color_RAGERS_x[i], color_RAGERS_y[i]])

# a = color_RAGERS_x
# b = color_RAGERS_y
# plt.scatter(a, b)

# c1 = color_RAGERS_err_x
# c2 = color_RAGERS_err_y[0]
 
# plt.errorbar(a, b, xerr=c1, fmt="o")
# plt.errorbar(a, b, yerr=c2, fmt="o")

# plt.show()




# plt.figure()

# a = color_RAGERS_x
# b = color_RAGERS_y
# plt.scatter(a, b)
 
# c = color_RAGERS_err_x
# d = color_RAGERS_err_y

# plt.errorbar(a, b, xerr=c, yerr=d, fmt="o", color="r")
 
# plt.show()

