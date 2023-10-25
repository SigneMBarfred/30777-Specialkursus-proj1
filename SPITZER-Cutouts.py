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
import matplotlib.colors as colors

#should this be a for loop?
# Open the FITS file and extract data and header
fits_file_path = 'C:/Users/signe/Documents/DTU/Specialkursus_radiogal/irac_map_bcd_pbcd/r10789376/ch1/pbcd/SPITZER_I1_10789376_0000_7_E8309859_maic.fits'
hdu_spitzer = fits.open(fits_file_path)
spitzer_data = hdu_spitzer[0].data
spitzer_hdr = hdu_spitzer[0].header
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
plt.title('Spitzer Image')
plt.xlabel('RA (J2000)')
plt.ylabel('DEC (J2000)')

# Read positions from your file
positions_file = 'C:/Users/signe/Documents/DTU/Specialkursus_radiogal/RAGERStxt.txt'
positions = np.loadtxt(positions_file, usecols=(1, 2), skiprows=1)  # Skip the header row
idno = np.loadtxt(positions_file,usecols=(0),skiprows=1)
# Define the cutout size in pixels

spitzer_pixelscale = u.pixel_scale(0.6*u.arcsec/u.pixel)

size = ((50*u.arcsec).to(u.pixel, spitzer_pixelscale), (50*u.arcsec).to(u.pixel, spitzer_pixelscale))  # the set size in arcsec is converted to pix
#2 buesec/arcsec pr. pixel is a true approx for SCUBA but is it usable here?


# Create contour levels to match the original plot
levels = (3, 4, 5, 6, 8, 10, 12, 14)

count=0

# Create a new figure for each cutout
for position in positions:
    loc = SkyCoord(position[0] * u.deg, position[1] * u.deg)
    cutout = Cutout2D(original_image, loc, size, wcs=wcs_spitzer)

    

    plt.figure()
    ax_cutout = plt.subplot(projection=cutout.wcs)
    im = ax_cutout.imshow(cutout.data, norm=norm, cmap='Greys')
    plt.colorbar(im, ax=ax_cutout)
    ax_cutout.grid(color='white', ls='solid')
    ax_cutout.set_xlabel('RA (J2000)')
    ax_cutout.set_ylabel('DEC (J2000)')
    ax_cutout.set_title('IRAC Ch1 map of Obj_id' + str(int(idno[count])) + ' \n at ' + str(position) + '\n')
    
    # Add the contour lines to the cutout plot

    # Write the cutout to a new FITS file
    cutout_filename = 'cutout_spitzer_ch1_id_' + str(int(idno[count])) + '.fits'
    fits.writeto(cutout_filename, cutout.data, overwrite=True)

    count += 1





