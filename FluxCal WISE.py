#%%

################
### Get Flux ###
################

import pandas as pd
import numpy as np
import csv
import math

from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
import astropy.units as u

# load Scuba data
# hdu_scuba = fits.open('C:/Users/Nikolaj Lange Dons/OneDrive - Danmarks Tekniske Universitet/Dokumenter/5 semester/Space special kursus/Saved files/cutout1.0.fits')
hdu_scuba = fits.open('C:/Users/olive/Dropbox/Group-1/Source Catalogs/maps/3C239_850_mf_cal_crop_snr.fits')
scuba850_data_cube = hdu_scuba[0].data  # This is a 3D cube
hdr_scuba = hdu_scuba[0].header
wcs_scuba = WCS(hdr_scuba)
df_scuba = pd.DataFrame(scuba850_data_cube[0,:,:]) # Pandas dataframe

# RAGERS koordinater txt fil
positions_file = 'C:/Users/olive/Dropbox/5. semester/Special kursus/RAGERStxt_uden_obj36.txt'
txt_coord = pd.read_csv(positions_file, sep=" ") # Pandas dataframe
obj_id = txt_coord['obj_id']

# WISE data
WISE = pd.read_csv('C:/Users/olive/Dropbox/5. semester/Special kursus/table_irsa_catalog_search_results.csv')
WISE_header = WISE.columns

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
    RAGERS[i] = pd.DataFrame(WISE[WISE.columns[0:36]][WISE[str(mask_id)]])
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

# Get only value from dataframe cell
# RAGERS[i]['designation'].values[0]

# SCUBA/RAGER
# z = 1.781

###################################
# Loop to calculate flux for WISE sources around RAGERS

for i in range(0,len(txt_coord['ra'])):
    m_vega_W1 = RAGERS[i]['w1mpro']
    m_vega_W2 = RAGERS[i]['w2mpro']
    m_vega_W3 = RAGERS[i]['w3mpro']
    m_vega_W4 = RAGERS[i]['w4mpro']

    # Flux uncerntanty in Vega Magnitude  units for different chanels
    m_vega_error_W1 = RAGERS[i]['w1sigmpro']
    m_vega_error_W2 = RAGERS[i]['w2sigmpro']
    m_vega_error_W3 = RAGERS[i]['w3sigmpro']
    m_vega_error_W4 = RAGERS[i]['w4sigmpro']

    # Constants
    Fv0_W1 = 309.540
    Fv0_W2 = 171.757
    Fv0_W3 = 31.678
    Fv0_W4 = 8.363

    # Equation 1 - Vega Magnitudes to Flux Density
    RAGERS[i]['flux_W1'] = Fv0_W1 * 10**(-m_vega_W1/2.5)
    RAGERS[i]['flux_W2'] = Fv0_W2 * 10**(-m_vega_W2/2.5)
    RAGERS[i]['flux_W3'] = Fv0_W3 * 10**(-m_vega_W3/2.5)
    RAGERS[i]['flux_W4'] = Fv0_W4 * 10**(-m_vega_W4/2.5)

    # Equation 1 - Vega Magnitudes uncertainty to Flux Density uncertainty
    RAGERS[i]['flux_error_W1'] = Fv0_W1 * 10**(-m_vega_error_W1/2.5)
    RAGERS[i]['flux_error_W2'] = Fv0_W2 * 10**(-m_vega_error_W2/2.5)
    RAGERS[i]['flux_error_W3'] = Fv0_W3 * 10**(-m_vega_error_W3/2.5)
    RAGERS[i]['flux_error_W4'] = Fv0_W4 * 10**(-m_vega_error_W4/2.5)

############################
### flux til sources     ###

for i in [0, 2, 5, 6, 8, 11, 12, 13, 14, 15, 17, 21]:
    flux_val1 = RAGERS[i]['flux_W1'].values[0]
    flux_val2 = RAGERS[i]['flux_W2'].values[0]
    flux_val3 = RAGERS[i]['flux_W3'].values[0]
    flux_val4 = RAGERS[i]['flux_W4'].values[0]
    print(f'{obj_id[i]}: Ch1 = {flux_val1}\n Ch2 = {flux_val2}\n Ch3 = {flux_val3}\n CH4 = {flux_val4}')

### RAGERS calling commands ###

# RAGERS[RAGER number with base zero][string column name]
# To call multiple columns:
    # RAGERS[15][['ang_dist', 'flux_W1']]

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

