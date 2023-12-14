# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 14:13:04 2023

@author: group1
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
import matplotlib.colors as colors

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
scuba_pixelscale = u.pixel_scale(2*u.arcsec/u.pixel)
(50*u.arcsec).to(u.pixel, scuba_pixelscale)
size = ((50*u.arcsec).to(u.pixel, scuba_pixelscale), (50*u.arcsec).to(u.pixel, scuba_pixelscale))  # the set size in arcsec is converted to pix


# Create contour levels (as defined in biggs 2010)
levels = (3, 4, 5, 6, 8, 10, 12, 14)

# Loop through positions and create cutouts
count = 0

for position in positions:
    
    #below is structured as in https://github.com/astropy/astropy/issues/6911 
    # Load the image and the world coordinates from the header
    hdu = fits.open('C:/Users/signe/Documents/DTU/Specialkursus_radiogal/3C239_850_mf_cal_crop_snr.fits')
    hduheader = hdu[0].header
    wcs = WCS(hduheader).dropaxis(2)
  
 
    # Making the cutout using the wcs
    loc = SkyCoord(position[0]*u.deg, position[1]*u.deg)
    cutout = Cutout2D(hdu[0].data[0], position=loc, size=size, wcs=wcs)
    
    # Update image data from the cutout
    hdu.data = cutout.data
    
    # Update the WCS from the cutout
    hdu[0].header.update(cutout.wcs.to_header())
    
    # Replace original image by the new "trimmed" version
    # Write the cutout to a new FITS file
    cutout_filename = 'scuba_cutout_'+str(idno[count])+'.fits'
    hdu.writeto(cutout_filename, overwrite=True)
    count += 1

    #plot the cutouts
    plt.figure()
    ax_cutout = plt.subplot(projection=cutout.wcs)
    im = ax_cutout.imshow(cutout.data, norm=norm, cmap='viridis')  # Use 'viridis' colormap
    plt.colorbar(im, ax=ax_cutout)  # Add colorbar to the current axes
    ax_cutout.grid(color='white', ls='solid')
    ax_cutout.set_xlabel('RA (J2000)')
    ax_cutout.set_ylabel('DEC (J2000)')
    ax_cutout.set_title('40\" by 40\" cutout @ 850 \u03BCm SCUBA map'+'\n Centered at '+str(position)+'\n')

    # Add the contour lines to the cutout plot
    ax_cutout.contour(cutout.data, levels=levels, colors='red', alpha=0.5)

# test that the generated fits files have correct wcs data
cutout_test = fits.open('C:/Users/signe/Documents/DTU/Specialkursus_radiogal/scuba_cutout_1.0.fits')
cutout_data = cutout_test[0].data
cutout_wcs = WCS(cutout_test[0].header)

#plot cutout example
plt.figure()
ax_cutout = plt.subplot(projection=cutout_wcs.dropaxis(2))
im = ax_cutout.imshow(cutout_data[0], norm=norm, cmap='viridis')
plt.colorbar(im, ax=ax_cutout)  # Add colorbar to the current axes
ax_cutout.grid(color='white', ls='solid')
ax_cutout.set_xlabel('RA (J2000)')
ax_cutout.set_ylabel('DEC (J2000)')
ax_cutout.set_title('test af cutout 1')
ax_cutout.contour(cutout_data[0], levels=levels, colors='red', alpha=0.5)

#issue w plot above: wcs naxis is 400 , 401 but should be 25,25(?) also does not correspond to the obj1 cutout plottet via forloop
