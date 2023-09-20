# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 12:32:44 2023

@author: signe
"""
#input: path to specific folder containing choosen channel of spitzer data 
#output: spitzer data stored in list

from astropy.io import fits
from astropy.wcs import WCS
import os,glob
from astropy.utils.data import get_pkg_data_filename


folder_path = 'C:/Users/signe/Documents/DTU/Specialkursus_radiogal/irac_map_bcd_pbcd/r10789376/ch2/bcd'


def spitzer_load(folder_path):
    n = 0 
    for filename in glob.glob(os.path.join(folder_path, '*.fits')):
        hdu = fits.open(filename)
        hdr = hdu[0].header  # the primary HDU header with metadata
        data = hdu[0].data
        #append coordinates to list? for all files in directory, so as to make crossref w ragers easier? 
        n = n+1
        wcs = WCS(hdr)
        
    print(n,'fits files found and read')
    
    if 'ch1' in folder_path:
        wavelength = 3.6*10**-6
        
    elif 'ch2' in folder_path:
        wavelength = 4.5*10**-6

    elif 'ch3' in folder_path:
         wavelength = 5.8*10**-6

    elif 'ch4' in folder_path:
        wavelength = 8.0*10**-6
      
    else: 
        wavelength = 'unknown'

    print('wavelength for measurements is', wavelength,'meters')

    return data,hdr,n,wcs

############ test ############

hdr = spitzer_load(folder_path)[1]
wcs = spitzer_load(folder_path)[3]

RA = hdr['RA_REF']
DEC = hdr['DEC_REF']

