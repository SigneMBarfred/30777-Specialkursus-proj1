# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 14:40:55 2023

@author: signe


A script that overlays the SCUBA 850micron snr map with the WISE ch1,2,3 and 4 map
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
#open wise cutout and scuba cutout


#first get the obj ids
positions_file = 'C:/Users/signe/Documents/DTU/Specialkursus_radiogal/RAGERStxt.txt'
positions = np.loadtxt(positions_file, usecols=(1, 2), skiprows=1)  # Skip the header row
idno = np.loadtxt(positions_file,usecols=(0),skiprows=1)

count = 0

for obj in idno:
    
    try:
        #obj is an id no.
        
        #wise file names are:
        #'cutout_wise_w[no]_id_' + obj + '.fits'
        hdu_wise = fits.open('C:/Users/signe/Documents/DTU/Specialkursus_radiogal/cutout_wise_w3_id_' 
                                + str(int(obj)) + '.fits')
        data_wise = hdu_wise[0].data
        hdr_wise = hdu_wise[0].header
        wcs_wise = WCS(hdr_wise)
        
        
        
        #scuba file names are:
        #'scuba_cutout_'+ obj +'.fits'
        hdu_scuba = fits.open('C:/Users/signe/Documents/DTU/Specialkursus_radiogal/scuba_cutout_'
                              + str(obj) +'.fits')
        scuba850_data_cube = hdu_scuba[0].data  # This is a 3D cube
        scuba_data = scuba850_data_cube
        hdr_scuba = hdu_scuba[0].header
        wcs_scuba = WCS(hdr_scuba)
        wcs_slice = wcs_scuba.dropaxis(2) 
    
        levels = (3, 4, 5, 6, 8, 10, 12, 14)
    
        loc = SkyCoord(positions[count][0] * u.deg, positions[count][1] * u.deg, frame=ICRS, obstime='J2000')
    
        beam = SphericalCircle(center=loc, radius = 18*u.arcsec, facecolor="None",  edgecolor="red", linewidth=2)
    
    
        vmin = np.nanmin(data_wise)
        vmax = np.nanmax(data_wise)
        norm = colors.SymLogNorm(linthresh=0.3, linscale=2.5, vmin=vmin, vmax=vmax)
        
        
        plt.figure()
        ax_wise = plt.subplot(projection=wcs_wise)
        im = ax_wise.imshow(data_wise, norm=norm, cmap='Greys')
        
        cb = plt.colorbar(im, ax=ax_wise)
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        cb.update_ticks()
        ax_wise.grid(color='black', ls='solid')
        ax_wise.set_xlabel('RA (J2000)')
        ax_wise.set_ylabel('DEC (J2000)')
        ra = ax_wise.coords[0]
        dec = ax_wise.coords[1]
        ra.set_ticks(number=4)
        dec.set_ticks(number=4)
        
        ax_wise.set_title('WISE W3 map of Obj_id' + str(int(obj)) )
        #ax_wise.contour(scuba_data, levels=levels, colors='blue',alpha=0.5)
      
        
        ax_wise.contour(scuba_data, transform=ax_wise.get_transform(wcs_slice),
                        levels=levels,
                        colors='blue', 
                        linewidths=0.5)
        #ax_wise.scatter(positions[count][0], positions[count][1], transform=ax_wise.get_transform('fk5'), s= 15*u.arcsec,
               #edgecolor='red', facecolor='none')
        
        
        #circle code based on https://gist.dropgithub.com/PaulHancock/7f7a4231d1ed9838775e965d90b25393 
        circle = SphericalCircle(loc, (15/3600)*u.deg, transform=ax_wise.get_transform('fk5'),color='red', alpha=0.2, facecolor='None')
        ax_wise.add_patch(circle)
        
        plt.savefig('scuba_overlay_wise_w3_id_' + str(int(idno[count])) + '.png')
        
    
    
    except FileNotFoundError:
        print("No cutout image of obj_id "+str(int(obj))+" as it lies outside WISE image area")
        pass
    
    count += 1

