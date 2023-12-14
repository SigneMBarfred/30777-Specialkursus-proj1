# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 14:45:47 2023

@author: olive & signe
"""

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
import matplotlib.colors as colors

#should this be a for loop?
# Open the FITS file and extract data and header
fits_file_path = 'C:/Users/signe/Documents/DTU/Specialkursus_radiogal/irac_map_bcd_pbcd/r10789376/ch3/pbcd/SPITZER_I3_10789376_0000_7_E8309809_maic.fits'
hdu_spitzer = fits.open(fits_file_path)
spitzer_data = hdu_spitzer[0].data
spitzer_hdr = hdu_spitzer[0].header
#print(spitzer_hdr)
wcs_spitzer = WCS(spitzer_hdr)



# Create a figure for the original Spitzer image
plt.figure()
vmin = np.nanmin(spitzer_data)
vmax = np.nanmax(spitzer_data)
norm = colors.SymLogNorm(linthresh=0.5, linscale=3, vmin=vmin, vmax=vmax)
original_image = spitzer_data
plt.subplot(projection=wcs_spitzer)
plt.imshow(original_image, origin='lower', norm=norm, cmap='viridis')
plt.colorbar()
plt.grid(color='white', ls='solid')
plt.title('Spitzer ch3 Image')
plt.xlabel('RA (J2000)')
plt.ylabel('DEC (J2000)')

# Read positions from your file
positions_file = 'C:/Users/signe/Documents/DTU/Specialkursus_radiogal/RAGERStxt.txt'
positions = np.loadtxt(positions_file, usecols=(1, 2), skiprows=1)  # Skip the header row
idno = np.loadtxt(positions_file,usecols=(0),skiprows=1)
# Define the cutout size in pixels

spitzer_pixelscale = u.pixel_scale(0.6*u.arcsec/u.pixel)

size = u.Quantity([50, 50], u.arcsec) # the set size in arcsec is converted to pix
#2 buesec/arcsec pr. pixel is a true approx for SCUBA but is it usable here?


# Create contour levels to match the original plot
levels = (3, 4, 5, 6, 8, 10, 12, 14)

count=0

# Create a new figure for each cutout
for position in positions[0:21]:
    loc = SkyCoord(position[0] * u.deg, position[1] * u.deg, frame=ICRS, obstime='J2000')
    cutout = Cutout2D(original_image, loc, size, wcs=wcs_spitzer)
    hdu_spitzer[0].data = cutout.data



    plt.figure()
    ax_cutout = plt.subplot(projection=cutout.wcs)
    im = ax_cutout.imshow(cutout.data, norm=norm, cmap='Greys')
    plt.colorbar(im, ax=ax_cutout)
    ax_cutout.grid(color='white', ls='solid')
    ax_cutout.set_xlabel('RA (J2000)')
    ax_cutout.set_ylabel('DEC (J2000)')
    ax_cutout.set_title('IRAC Ch3 map of Obj_id' + str(int(idno[count])) + ' \n at ' + str(position) + '\n')
    
    plt.savefig('cutout_spitzer_ch3_id_' + str(int(idno[count])) + '.png')
    
    #update header
    hdu_spitzer[0].header.update(cutout.wcs.to_header()) 

    # Write the cutout to a new FITS file
    cutout_filename = 'cutout_spitzer_ch3_id_' + str(int(idno[count])) + '.fits'
    fits.writeto(cutout_filename, cutout.data, header=hdu_spitzer[0].header, overwrite=True)

    
    count += 1



# test that the generated fits files have correct wcs data
cutout_test = fits.open('C:/Users/signe/Documents/DTU/Specialkursus_radiogal/cutout_spitzer_ch3_id_1.fits')
cutout_data = cutout_test[0].data
cutout_wcs = WCS(cutout_test[0].header)

#plot cutout example
plt.figure()
ax_cutout = plt.subplot(projection=cutout_wcs)
im = ax_cutout.imshow(cutout_data, norm=norm, cmap='viridis')
plt.colorbar(im, ax=ax_cutout)  # Add colorbar to the current axes
ax_cutout.grid(color='white', ls='solid')
ax_cutout.set_xlabel('RA (J2000)')
ax_cutout.set_ylabel('DEC (J2000)')
ax_cutout.set_title('test af spitzer ch3 cutout 1')
ax_cutout.contour(cutout_data, levels=levels, colors='red', alpha=0.5)


