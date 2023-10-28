# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 12:06:20 2023

@author: Nikolaj Lange Dons
"""

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
hdu_scuba = fits.open('C:/Users/Nikolaj Lange Dons/OneDrive - Danmarks Tekniske Universitet/Dokumenter/5 semester/Space special kursus/Saved files/cutout1.0.fits')
scuba850_data_cube = hdu_scuba[0].data  # This is a 3D cube
hdr_scuba = hdu_scuba[0].header
wcs_scuba = WCS(hdr_scuba)
df_scuba = pd.DataFrame(scuba850_data_cube[0,:,:]) # Pandas dataframe

# RAGERS koordinater txt fil
positions_file = 'C:/Users/Nikolaj Lange Dons/OneDrive - Danmarks Tekniske Universitet/Dokumenter/5 semester/Space special kursus/RAGERStxt.txt'
txt_coord = pd.read_csv(positions_file, sep=" ") # Pandas dataframe
obj_id = txt_coord['obj_id']
# np_coord = pd.DataFrame(txt_coord).to_numpy()
# RAGERS_coord = np.transpose(np.array([np_coord[:,1],np_coord[:,2]]))

# WISE data
WISE = pd.read_csv('C:/Users/Nikolaj Lange Dons/OneDrive - Danmarks Tekniske Universitet/Dokumenter/5 semester/Space special kursus/Saved files/table_irsa_catalog_search_results.csv')
WISE_header = WISE.columns
# np_WISE = pd.DataFrame(WISE).to_numpy()
# ra = np_WISE[:,1] 
# dec = np_WISE[:,2]
# wise = np.asarray(np_WISE[:,1:], dtype = np.float64)


####################################
## Distance to RAGERS

# def distance_to_centrum(df, center_x, center_y):
#     df['distance_to_centrum'] = np.sqrt((df['ra'] - center_x)**2 + (df['dec'] - center_y)**2)
#     return df

# x_coord = txt_coord['ra']
# y_coord = txt_coord['dec']
# WISE_coord = WISE[['ra','dec']]

# WISE_mask = distance_to_centrum(WISE_coord, x_coord[0], y_coord[0])

#####################
# Find target sources around RAGERS

# Mask based on distance to RAGERS
def distance_to_centrum(df, center_x, center_y, radius):
    df['distance_to_centrum'] = np.sqrt((df['ra'] - center_x)**2 + (df['dec'] - center_y)**2)
    df['mask'] = df['distance_to_centrum'] <= radius
    return df['mask']

# Define coordinates
x_coord = txt_coord['ra']
y_coord = txt_coord['dec']
WISE_coord = WISE[['ra','dec']]

# Loop to ceate new columns with masks from function
for i in range(0,len(txt_coord['ra'])):
    mask_name = 'mask'
    mask_id = (mask_name + str(i+1))
    WISE[str(mask_name + str(i+1))] = distance_to_centrum(WISE_coord, x_coord[i], y_coord[i], 15.0/3600.0)

# n dataframe navne
RAGERS = [['RAGER'+ str(_dummy)] for _dummy in range(1, 25)]

# Loop to create list of dataframes with sources within 15 arcsec of RAGERS
for i in range(0,len(txt_coord['ra'])):
    mask_name = 'mask'
    mask_id = (mask_name + str(i+1))
    RAGERS[i] = pd.DataFrame(WISE[WISE.columns[0:36]][WISE[str(mask_name + str(i+1))]])

# Count sources around each RAGER
n_source = np.array(np.zeros(len(txt_coord['ra'])))
for i in range(0,len(txt_coord['ra'])):
    n_source[i] = np.array([len(RAGERS[i]['w1mpro'])])


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


##############################################
### Calculate flux for all of WISE sources ###

# Vega Magnitudes for different chanels
m_vega_W1 = WISE['w1mpro']
m_vega_W2 = WISE['w2mpro']
m_vega_W3 = WISE['w3mpro']
m_vega_W4 = WISE['w4mpro']

# Flux uncerntanty in Vega Magnitude units for different chanels
m_vega_error_W1 = WISE['w1sigmpro']
m_vega_error_W2 = WISE['w2sigmpro']
m_vega_error_W3 = WISE['w3sigmpro']
m_vega_error_W4 = WISE['w4sigmpro']

# Constants
Fv0_W1 = 309.540
Fv0_W2 = 171.757
Fv0_W3 = 31.678
Fv0_W4 = 8.363

# Equation 1 - Vega Magnitudes to Flux Density
flux_W1 = Fv0_W1 * 10**(-m_vega_W1/2.5)
flux_W2 = Fv0_W2 * 10**(-m_vega_W2/2.5)
flux_W3 = Fv0_W3 * 10**(-m_vega_W3/2.5)
flux_W4 = Fv0_W4 * 10**(-m_vega_W4/2.5)

# Equation 1 - Vega Magnitudes uncertainty to Flux Density uncertainty
flux_error_W1 = Fv0_W1 * 10**(-m_vega_error_W1/2.5)
flux_error_W2 = Fv0_W2 * 10**(-m_vega_error_W2/2.5)
flux_error_W3 = Fv0_W3 * 10**(-m_vega_error_W3/2.5)
flux_error_W4 = Fv0_W4 * 10**(-m_vega_error_W4/2.5)



