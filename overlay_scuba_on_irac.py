# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 14:09:38 2023

@author: signe

A script that overlays the SCUBA 850micron snr map with the IRAC ch1,2,3 and 4 map
Aim: Visualize possible SMG counterparts
"""

from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
# for each obj id 
#open spitzer cutout and scuba cutout


#first get the obj ids
positions_file = 'C:/Users/signe/Documents/DTU/Specialkursus_radiogal/RAGERStxt.txt'
positions = np.loadtxt(positions_file, usecols=(1, 2), skiprows=1)  # Skip the header row
idno = np.loadtxt(positions_file,usecols=(0),skiprows=1)

count = 0

for obj in idno:
    #obj is an id no.
    
    #irac file names are:
    #'cutout_spitzer_ch1_id_' + obj + '.fits'
    hdu_spitzer = fits.open('C:/Users/signe/Documents/DTU/Specialkursus_radiogal/cutout_spitzer_ch1_id_' + str(int(obj)) + '.fits')
    data_spitzer = hdu_spitzer[0].data
    hdr_spitzer = hdu_spitzer[0].header
    wcs_spitzer = WCS(hdr_spitzer)
    
    
    
    #scuba file names are:
    #'scuba_cutout_'+ obj +'.fits'
    hdu_scuba = fits.open('C:/Users/signe/Documents/DTU/Specialkursus_radiogal/scuba_cutout_'+ str(obj) +'.fits')
    scuba850_data_cube = hdu_scuba[0].data  # This is a 3D cube
    scuba_data = scuba850_data_cube[0]
    hdr_scuba = hdu_scuba[0].header
    wcs_scuba = WCS(hdr_scuba)
    wcs_slice = wcs_scuba.dropaxis(2) 

    levels = (3, 4, 5, 6, 8, 10, 12, 14)


    vmin = np.nanmin(data_spitzer)
    vmax = np.nanmax(data_spitzer)
    norm = colors.SymLogNorm(linthresh=0.5, linscale=3, vmin=vmin, vmax=vmax)
    plt.figure()
    ax_spitzer = plt.subplot(projection=wcs_spitzer)
    im = ax_spitzer.imshow(data_spitzer, norm=norm, cmap='Greys')
    plt.colorbar(im, ax=ax_spitzer)
    ax_spitzer.grid(color='white', ls='solid')
    ax_spitzer.set_xlabel('RA (J2000)')
    ax_spitzer.set_ylabel('DEC (J2000)')
    ax_spitzer.set_title('IRAC Ch1 map of Obj_id' + str(int(obj)) + ' \n at ' + str(positions[count]) + '\n')
    ax_spitzer.contour(scuba_data,transform=ax_spitzer.get_transform(wcs_slice), levels=levels, colors='blue',alpha=0.5)
    #ax.contour(hdu_scuba.data, vmin=-1, vmax=10, transform=ax.get_transform(WCS(hdu_scuba.header)),
           #levels=[3, 4,5,6,7,8, 9,10,11,12], colors='White', linewidths=0.5)
    
    count += 1
#show spitzer cutout
# overlay contours from SCUBA
