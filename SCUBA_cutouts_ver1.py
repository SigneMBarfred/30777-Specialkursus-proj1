# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 14:10:04 2023

@author: group1
IT WORKS
code that reads scuba data, plots entire picture and thereafter creates cutouts at SMG-positions 
and saves those as new fits files

"""

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
import matplotlib.colors as colors
from astropy.visualization.wcsaxes import SphericalCircle
from matplotlib.collections import PatchCollection 

# Open the FITS file and extract data and header
hdu_scuba = fits.open('C:/Users/signe/Documents/DTU/Specialkursus_radiogal/3C239_850_mf_cal_crop_snr.fits')
scuba850_data_cube = hdu_scuba[0].data  # This is a 3D cube
hdr_scuba = hdu_scuba[0].header
wcs_scuba = WCS(hdr_scuba)

# Create a figure for the original image
plt.figure()
slice_index = 0  # Change this to the desired slice index within the cube
wcs_slice = wcs_scuba.dropaxis(2)  # Remove the third dimension for visualization
vmin = np.nanmin(scuba850_data_cube[slice_index])
vmax = np.nanmax(scuba850_data_cube[slice_index])
norm = colors.SymLogNorm(linthresh=3, linscale=100, vmin=vmin, vmax=vmax)
original_image = scuba850_data_cube[slice_index]


#plot entire image
plt.subplot(projection=wcs_slice)
plt.imshow(original_image, origin='lower', norm=norm, cmap='viridis')  # Use 'viridis' colormap
plt.colorbar()
plt.grid(color='white', ls='solid')
plt.title('850' + r'$\mu$' + 'm')
plt.xlabel('RA (J2000)')
plt.ylabel('DEC (J2000)')

# Read positions from a file where x-coordinates are in the second column and y-coordinates are in the third column
positions_file = 'C:/Users/signe/Documents/DTU/Specialkursus_radiogal/RAGERStxt.txt'
positions = np.loadtxt(positions_file, usecols=(1, 2), skiprows=1)  # Skip the header row
idno = np.loadtxt(positions_file,usecols=(0),skiprows=1)

# Define the cutout size in pixels
#2 arcsec pr. pixel
size = u.Quantity([50, 50], u.arcsec)
# Create contour levels (as defined in biggs 2010)
levels = (3, 4, 5, 6, 8, 10, 12, 14)

# Loop through positions and create cutouts
count = 0 
    
for position in positions:
    loc = SkyCoord(position[0]*u.deg, position[1]*u.deg, frame=FK5, obstime='J2000')
    cutout = Cutout2D(original_image, loc, size, wcs=wcs_slice)
    #cutout.writeto('cutout'+str(idno[count]).fits) #writes a fits file for each cutout scuba data
    
    
    # Put the cutout image in the FITS HDU
    hdu_scuba.data = cutout.data

    # Update the FITS header with the cutout WCS
    new_hdr = hdr_scuba.copy
    
    hdu_scuba[0].header.update(cutout.wcs.to_header()) 

    # Write the cutout to a new FITS file
    cutout_filename = 'scuba_cutout_'+str(idno[count])+'.fits'
    fits.writeto(cutout_filename, cutout.data, header=hdu_scuba[0].header, overwrite=True)

    count += 1
    
    
    
    
    # Create a new figure for each cutout
    plt.figure()
    ax_cutout = plt.subplot(projection=cutout.wcs)
    im = ax_cutout.imshow(cutout.data, norm=norm, cmap='viridis')  # Use 'viridis' colormap
    plt.colorbar(im, ax=ax_cutout)  # Add colorbar to the current axes
    ax_cutout.grid(color='white', ls='solid')
    ax_cutout.set_xlabel('RA (J2000)')
    ax_cutout.set_ylabel('DEC (J2000)')
    ax_cutout.set_title('50\" by 50\" cutout @ 850 \u03BCm SCUBA map'+'\n Centered at '+str(position)+'\n')

    # Add the contour lines to the cutout plot
    ax_cutout.contour(cutout.data, levels=levels, colors='red', alpha=0.5)

plt.show()

# test that the generated fits files have correct wcs data
cutout_test = fits.open('C:/Users/signe/Documents/DTU/Specialkursus_radiogal/scuba_cutout_1.0.fits')
cutout_data = cutout_test[0].data
cutout_wcs = WCS(cutout_test[0].header)
cutout_wcs = cutout_wcs.dropaxis(2)
#plot cutout example
plt.figure()
ax_cutout = plt.subplot(projection=cutout_wcs)
im = ax_cutout.imshow(cutout_data, norm=norm, cmap='viridis')
plt.colorbar(im, ax=ax_cutout)  # Add colorbar to the current axes
ax_cutout.grid(color='white', ls='solid')
ax_cutout.set_xlabel('RA (J2000)')
ax_cutout.set_ylabel('DEC (J2000)')
ax_cutout.set_title('test af cutout 1 v1')
ax_cutout.contour(cutout_data, levels=levels, colors='red', alpha=0.5)
