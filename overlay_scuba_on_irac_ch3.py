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
from matplotlib.patches import Circle 
from astropy.visualization.wcsaxes import SphericalCircle
from matplotlib.collections import PatchCollection 
from matplotlib import cm 
from astropy.coordinates import ICRS, Galactic, FK4, FK5
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
# for each obj id 
#open spitzer cutout and scuba cutout


#first get the obj ids
positions_file = 'C:/Users/signe/Documents/DTU/Specialkursus_radiogal/RAGERStxt.txt'
positions = np.loadtxt(positions_file, usecols=(1, 2), skiprows=1)  # Skip the header row
idno = np.loadtxt(positions_file,usecols=(0),skiprows=1)

count = 0

for obj in idno:
    
    try:
        #obj is an id no.
        
        #irac file names are:
        #'cutout_spitzer_ch1_id_' + obj + '.fits'
        hdu_spitzer = fits.open('C:/Users/signe/Documents/DTU/Specialkursus_radiogal/cutout_spitzer_ch3_id_' 
                                + str(int(obj)) + '.fits')
        data_spitzer = hdu_spitzer[0].data
        hdr_spitzer = hdu_spitzer[0].header
        wcs_spitzer = WCS(hdr_spitzer)
        
        
        
        #scuba file names are:
        #'scuba_cutout_'+ obj +'.fits'
        hdu_scuba = fits.open('C:/Users/signe/Documents/DTU/Specialkursus_radiogal/scuba_cutout_'
                              + str(obj) +'.fits')
        scuba850_data_cube = hdu_scuba[0].data  # This is a 3D cube
        scuba_data = scuba850_data_cube
        # scuba_data = scuba_data[0,:,:]
        hdr_scuba = hdu_scuba[0].header
        wcs_scuba = WCS(hdr_scuba)
        wcs_slice = wcs_scuba.dropaxis(2) 
    
        levels = np.array([3, 4, 5, 6, 8, 10, 12, 14])
    

    
    
        loc = SkyCoord(positions[count][0] * u.deg, positions[count][1] * u.deg, frame=ICRS, obstime='J2000')
    
        beam = SphericalCircle(center=loc, radius = 15*u.arcsec, facecolor="None",  edgecolor="red", linewidth=2)
    
    
        vmin = np.nanmin(data_spitzer)
        vmax = np.nanmax(data_spitzer)
        norm = colors.SymLogNorm(linthresh=1, linscale=2.5, vmin=vmin, vmax=vmax)
        
        
        
        plt.figure()
        ax_spitzer = plt.subplot(projection=wcs_spitzer)
        im = ax_spitzer.imshow(data_spitzer, norm=norm, cmap='Greys')
        
        cb = plt.colorbar(im, ax=ax_spitzer)
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        cb.update_ticks()
        ax_spitzer.grid(color='black', ls='solid')
        ax_spitzer.set_xlabel('RA (J2000)')
        ax_spitzer.set_ylabel('DEC (J2000)')
        ra = ax_spitzer.coords[0]
        dec = ax_spitzer.coords[1]
        plt.tick_params(axis='both', nbins=4) 
        ax_spitzer.coords[0].set_ticks(spacing=1. * u.arcsec)
        ra.set_ticks(number=4)
        dec.set_ticks(number=4)
        
        ax_spitzer.set_title('IRAC Ch3 map of Obj_id' + str(int(obj)))
        #ax_spitzer.contour(scuba_data, levels=levels, colors='blue',alpha=0.5)
      
        
        ax_spitzer.contour(scuba_data,
                           transform=ax_spitzer.get_transform(wcs_slice),
                           levels=levels,
                           colors='blue', 
                           linewidths=0.5)
        # ax_spitzer.scatter(positions[count][0], positions[count][1], transform=ax_spitzer.get_transform('fk5'), s= 18*u.arcsec,
        #        edgecolor='red', facecolor='none')
        
        circle = SphericalCircle(loc, (15/3600)*u.deg, transform=ax_spitzer.get_transform('fk5'),color='red', alpha=0.2, facecolor='None')
        ax_spitzer.add_patch(circle)
        plt.savefig('scuba_overlay_spitzer_ch3_id_' + str(int(idno[count])) + '.png')
        
    except FileNotFoundError:
        print("No cutout image of obj_id "+str(int(obj))+" as it lies outside IRAC image area")
        pass
    
    count += 1

