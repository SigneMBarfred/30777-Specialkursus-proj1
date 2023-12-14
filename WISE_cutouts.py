# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 14:04:40 2023
Creates WISE cutouts around SMGs
WISE has 4 bands and thus 4*23 cutouts are created
@author: signe
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



"""
W1

First part creates cutouts for W1 (3.4 microns) (mid-infrared)
Could probably easily be made into loop for all 4 passbands 

update variables fits_file_path and positions_file to run
"""

# Open the FITS file and extract data and header
fits_file_path = 'C:/Users/signe/Documents/DTU/Specialkursus_radiogal/WISE_all_bands_fits/1527p469_ac51-w1-int-3_ra152.9392302_dec46.47215342_asec600.000.fits'
hdu_wise = fits.open(fits_file_path)
wise_data = hdu_wise[0].data
wise_hdr = hdu_wise[0].header

wcs_wise = WCS(wise_hdr)



# Create a figure for the original Spitzer image
plt.figure()
vmin = np.nanmin(wise_data)
vmax = np.nanmax(wise_data)
norm = colors.SymLogNorm(linthresh=8, linscale=10, vmin=vmin, vmax=vmax)
original_image = wise_data
plt.subplot(projection=wcs_wise)
plt.imshow(original_image, origin='lower', norm=norm, cmap='viridis')
plt.colorbar()
plt.grid(color='white', ls='solid')
plt.title('WISE W1 Image of cluster 3C239 \n 600\u2033*600\u2033 at 3.4$\mu$m')
plt.xlabel('RA (J2000)')
plt.ylabel('DEC (J2000)')

# Read positions from your file
positions_file = 'C:/Users/signe/Documents/DTU/Specialkursus_radiogal/RAGERStxt.txt'
positions = np.loadtxt(positions_file, usecols=(1, 2), skiprows=1)  # Skip the header row
idno = np.loadtxt(positions_file,usecols=(0),skiprows=1)
# Define the cutout size in pixels

size = u.Quantity([50, 50], u.arcsec)
#2 buesec/arcsec pr. pixel is a true approx for SCUBA but is it usable here?


# Create contour levels to match the original plot
levels = (3, 4, 5, 6, 8, 10, 12, 14)

count=0

# Create a new figure for each cutout
for position in positions:
    
    try:
        loc = SkyCoord(position[0] * u.deg, position[1] * u.deg, frame=ICRS, obstime='J2000')
        cutout = Cutout2D(original_image, loc, size, wcs=wcs_wise, mode="partial", fill_value=np.nan)
        hdu_wise[0].data = cutout.data
    
    
    
        # plt.figure()
        # ax_cutout = plt.subplot(projection=cutout.wcs)
        # im = ax_cutout.imshow(cutout.data, norm=norm, cmap='Greys')
        # plt.colorbar(im, ax=ax_cutout)
        # ax_cutout.grid(color='black', ls='solid')
        # ax_cutout.set_xlabel('RA (J2000)')
        # ax_cutout.set_ylabel('DEC (J2000)')
        # ax_cutout.set_title('WISE Ch1 map of Obj_id' + str(int(idno[count])) + ' \n at ' + str(position) + '\n')
        
        
        #update header
        hdu_wise[0].header.update(cutout.wcs.to_header()) 
    
        # Write the cutout to a new FITS file
        cutout_filename = 'cutout_wise_w1_id_' + str(int(idno[count])) + '.fits'
        fits.writeto(cutout_filename, cutout.data, header=hdu_wise[0].header, overwrite=True)

    except:
        print("obj_id "+str(int(idno[count]))+" is outside WISE W1 image.")
        pass
    
    count += 1



# # test that the generated fits files have correct wcs data
# cutout_test = fits.open('C:/Users/signe/Documents/DTU/Specialkursus_radiogal/cutout_spitzer_ch1_id_1.fits')
# cutout_data = cutout_test[0].data
# cutout_wcs = WCS(cutout_test[0].header)

# #plot cutout example
# plt.figure()
# ax_cutout = plt.subplot(projection=cutout_wcs)
# im = ax_cutout.imshow(cutout_data, norm=norm, cmap='viridis')
# plt.colorbar(im, ax=ax_cutout)  # Add colorbar to the current axes
# ax_cutout.grid(color='white', ls='solid')
# ax_cutout.set_xlabel('RA (J2000)')
# ax_cutout.set_ylabel('DEC (J2000)')
# ax_cutout.set_title('test af wise cutout 1')
# ax_cutout.contour(cutout_data, levels=levels, colors='red', alpha=0.5)


"""
W2

Second part (below) creates cutouts for W2 (4.6 microns) (mid-far-infrared)
Could probably easily be made into loop for all 4 passbands 

update variables fits_file_path and positions_file to run
"""

# Open the FITS file and extract data and header
fits_file_path = 'C:/Users/signe/Documents/DTU/Specialkursus_radiogal/WISE_all_bands_fits/1527p469_ac51-w2-int-3_ra152.9392302_dec46.47215342_asec600.000.fits'
hdu_wise = fits.open(fits_file_path)
wise_data = hdu_wise[0].data
wise_hdr = hdu_wise[0].header

wcs_wise = WCS(wise_hdr)



# Create a figure for the original wise image
plt.figure()
vmin = np.nanmin(wise_data)
vmax = np.nanmax(wise_data)
norm = colors.SymLogNorm(linthresh=8, linscale=20, vmin=vmin, vmax=vmax)
original_image = wise_data
plt.subplot(projection=wcs_wise)
plt.imshow(original_image, origin='lower', norm=norm, cmap='viridis')
plt.colorbar()
plt.grid(color='white', ls='solid')
plt.title('WISE W2 Image of cluster 3C239 \n 600\u2033*600\u2033 at 4.6$\mu$m')
plt.xlabel('RA (J2000)')
plt.ylabel('DEC (J2000)')

# Read positions from your file
positions_file = 'C:/Users/signe/Documents/DTU/Specialkursus_radiogal/RAGERStxt.txt'
positions = np.loadtxt(positions_file, usecols=(1, 2), skiprows=1)  # Skip the header row
idno = np.loadtxt(positions_file,usecols=(0),skiprows=1)
# Define the cutout size in pixels

size = u.Quantity([50, 50], u.arcsec)
#2 buesec/arcsec pr. pixel is a true approx for SCUBA but is it usable here?


# Create contour levels to match the original plot
levels = (3, 4, 5, 6, 8, 10, 12, 14)

count=0

# Create a new figure for each cutout
for position in positions:
    
    try:
        loc = SkyCoord(position[0] * u.deg, position[1] * u.deg, frame=ICRS, obstime='J2000')
        cutout = Cutout2D(original_image, loc, size, wcs=wcs_wise, mode="partial", fill_value=np.nan)
        hdu_wise[0].data = cutout.data
    
    
    
        # plt.figure()
        # ax_cutout = plt.subplot(projection=cutout.wcs)
        # im = ax_cutout.imshow(cutout.data, norm=norm, cmap='Greys')
        # plt.colorbar(im, ax=ax_cutout)
        # ax_cutout.grid(color='black', ls='solid')
        # ax_cutout.set_xlabel('RA (J2000)')
        # ax_cutout.set_ylabel('DEC (J2000)')
        # ax_cutout.set_title('WISE W2 map of Obj_id' + str(int(idno[count])) + ' \n at ' + str(position) + '\n')
        
        
        #update header
        hdu_wise[0].header.update(cutout.wcs.to_header()) 
    
        # Write the cutout to a new FITS file
        cutout_filename = 'cutout_wise_w2_id_' + str(int(idno[count])) + '.fits'
        fits.writeto(cutout_filename, cutout.data, header=hdu_wise[0].header, overwrite=True)

    except:
        print("obj_id "+str(int(idno[count]))+" is outside WISE W2 image.")
        pass
    
    count += 1



# # test that the generated fits files have correct wcs data
# cutout_test = fits.open('C:/Users/signe/Documents/DTU/Specialkursus_radiogal/cutout_wise_w2_id_1.fits')
# cutout_data = cutout_test[0].data
# cutout_wcs = WCS(cutout_test[0].header)

# #plot cutout example
# plt.figure()
# ax_cutout = plt.subplot(projection=cutout_wcs)
# im = ax_cutout.imshow(cutout_data, norm=norm, cmap='viridis')
# plt.colorbar(im, ax=ax_cutout)  # Add colorbar to the current axes
# ax_cutout.grid(color='white', ls='solid')
# ax_cutout.set_xlabel('RA (J2000)')
# ax_cutout.set_ylabel('DEC (J2000)')
# ax_cutout.set_title('test af wise w2 cutout 1')
# ax_cutout.contour(cutout_data, levels=levels, colors='red', alpha=0.5)


"""
W3

Third part (below) creates cutouts for W3 (11.6 microns) (far-infrared)
Could probably easily be made into loop for all 4 passbands 

update variables fits_file_path and positions_file to run
"""

# Open the FITS file and extract data and header
fits_file_path = 'C:/Users/signe/Documents/DTU/Specialkursus_radiogal/WISE_all_bands_fits/1527p469_ac51-w3-int-3_ra152.9392302_dec46.47215342_asec600.000.fits'
hdu_wise = fits.open(fits_file_path)
wise_data = hdu_wise[0].data
wise_hdr = hdu_wise[0].header

wcs_wise = WCS(wise_hdr)



# Create a figure for the original wise image
plt.figure()
vmin = np.nanmin(wise_data)
vmax = np.nanmax(wise_data)
norm = colors.SymLogNorm(linthresh=0.001, linscale=0.01, vmin=vmin, vmax=vmax)
original_image = wise_data
plt.subplot(projection=wcs_wise)
plt.imshow(original_image, origin='lower', norm=norm, cmap='viridis')
plt.colorbar()
plt.grid(color='white', ls='solid')
plt.title('WISE W3 Image of cluster 3C239 \n 600\u2033*600\u2033 at 11.6$\mu$m')
plt.xlabel('RA (J2000)')
plt.ylabel('DEC (J2000)')

# Read positions from your file
positions_file = 'C:/Users/signe/Documents/DTU/Specialkursus_radiogal/RAGERStxt.txt'
positions = np.loadtxt(positions_file, usecols=(1, 2), skiprows=1)  # Skip the header row
idno = np.loadtxt(positions_file,usecols=(0),skiprows=1)
# Define the cutout size in pixels

size = u.Quantity([50, 50], u.arcsec)
#2 buesec/arcsec pr. pixel is a true approx for SCUBA but is it usable here?


# Create contour levels to match the original plot
levels = (3, 4, 5, 6, 8, 10, 12, 14)

count=0

# Create a new figure for each cutout
for position in positions:
    
    try:
        loc = SkyCoord(position[0] * u.deg, position[1] * u.deg, frame=ICRS, obstime='J2000')
        cutout = Cutout2D(original_image, loc, size, wcs=wcs_wise, mode="partial", fill_value=np.nan)
        hdu_wise[0].data = cutout.data
    
    
    
        # plt.figure()
        # ax_cutout = plt.subplot(projection=cutout.wcs)
        # im = ax_cutout.imshow(cutout.data, norm=norm, cmap='Greys')
        # plt.colorbar(im, ax=ax_cutout)
        # ax_cutout.grid(color='black', ls='solid')
        # ax_cutout.set_xlabel('RA (J2000)')
        # ax_cutout.set_ylabel('DEC (J2000)')
        # ax_cutout.set_title('WISE W3 map of Obj_id' + str(int(idno[count])) + ' \n at ' + str(position) + '\n')
        
        
        #update header
        hdu_wise[0].header.update(cutout.wcs.to_header()) 
    
        # Write the cutout to a new FITS file
        cutout_filename = 'cutout_wise_w3_id_' + str(int(idno[count])) + '.fits'
        fits.writeto(cutout_filename, cutout.data, header=hdu_wise[0].header, overwrite=True)

    except:
        print("obj_id "+str(int(idno[count]))+" is outside WISE W3 image.")
        pass
    
    count += 1



# # test that the generated fits files have correct wcs data
# cutout_test = fits.open('C:/Users/signe/Documents/DTU/Specialkursus_radiogal/cutout_wise_w3_id_1.fits')
# cutout_data = cutout_test[0].data
# cutout_wcs = WCS(cutout_test[0].header)

# #plot cutout example
# plt.figure()
# ax_cutout = plt.subplot(projection=cutout_wcs)
# im = ax_cutout.imshow(cutout_data, norm=norm, cmap='viridis')
# plt.colorbar(im, ax=ax_cutout)  # Add colorbar to the current axes
# ax_cutout.grid(color='white', ls='solid')
# ax_cutout.set_xlabel('RA (J2000)')
# ax_cutout.set_ylabel('DEC (J2000)')
# ax_cutout.set_title('test af wise w3 cutout 1')
# ax_cutout.contour(cutout_data, levels=levels, colors='red', alpha=0.5)



"""
W4

 Fourth part (below) creates cutouts for W4 (22 microns) (far-infrared)
Could probably easily be made into loop for all 4 passbands 

update variables fits_file_path and positions_file to run
"""

# Open the FITS file and extract data and header
fits_file_path = 'C:/Users/signe/Documents/DTU/Specialkursus_radiogal/WISE_all_bands_fits/1527p469_ac51-w4-int-3_ra152.9392302_dec46.47215342_asec600.000.fits'
hdu_wise = fits.open(fits_file_path)
wise_data = hdu_wise[0].data
wise_hdr = hdu_wise[0].header

wcs_wise = WCS(wise_hdr)



# Create a figure for the original Spitzer image
plt.figure()
vmin = np.nanmin(wise_data)
vmax = np.nanmax(wise_data)
norm = colors.SymLogNorm(linthresh=8, linscale=15, vmin=vmin, vmax=vmax)
original_image = wise_data
plt.subplot(projection=wcs_wise)
plt.imshow(original_image, origin='lower', norm=norm, cmap='viridis')
plt.colorbar()
plt.grid(color='white', ls='solid')
plt.title('WISE W4 Image of cluster 3C239 \n 600\u2033*600\u2033 at 22$\mu$m')
plt.xlabel('RA (J2000)')
plt.ylabel('DEC (J2000)')

# Read positions from your file
positions_file = 'C:/Users/signe/Documents/DTU/Specialkursus_radiogal/RAGERStxt.txt'
positions = np.loadtxt(positions_file, usecols=(1, 2), skiprows=1)  # Skip the header row
idno = np.loadtxt(positions_file,usecols=(0),skiprows=1)
# Define the cutout size in pixels

size = u.Quantity([50, 50], u.arcsec)
#2 buesec/arcsec pr. pixel is a true approx for SCUBA but is it usable here?


# Create contour levels to match the original plot
levels = (3, 4, 5, 6, 8, 10, 12, 14)

count=0

# Create a new figure for each cutout
for position in positions:
    
    try:
        loc = SkyCoord(position[0] * u.deg, position[1] * u.deg, frame=ICRS, obstime='J2000')
        cutout = Cutout2D(original_image, loc, size, wcs=wcs_wise, mode="partial", fill_value=np.nan)
        hdu_wise[0].data = cutout.data
    
    
    
        # plt.figure()
        # ax_cutout = plt.subplot(projection=cutout.wcs)
        # im = ax_cutout.imshow(cutout.data, norm=norm, cmap='Greys')
        # plt.colorbar(im, ax=ax_cutout)
        # ax_cutout.grid(color='black', ls='solid')
        # ax_cutout.set_xlabel('RA (J2000)')
        # ax_cutout.set_ylabel('DEC (J2000)')
        # ax_cutout.set_title('WISE W3 map of Obj_id' + str(int(idno[count])) + ' \n at ' + str(position) + '\n')
        
        
        #update header
        hdu_wise[0].header.update(cutout.wcs.to_header()) 
    
        # Write the cutout to a new FITS file
        cutout_filename = 'cutout_wise_w4_id_' + str(int(idno[count])) + '.fits'
        fits.writeto(cutout_filename, cutout.data, header=hdu_wise[0].header, overwrite=True)

    except:
        print("obj_id "+str(int(idno[count]))+" is outside WISE W4 image.")
        pass
    
    count += 1



# # test that the generated fits files have correct wcs data
# cutout_test = fits.open('C:/Users/signe/Documents/DTU/Specialkursus_radiogal/cutout_wise_w3_id_1.fits')
# cutout_data = cutout_test[0].data
# cutout_wcs = WCS(cutout_test[0].header)

# #plot cutout example
# plt.figure()
# ax_cutout = plt.subplot(projection=cutout_wcs)
# im = ax_cutout.imshow(cutout_data, norm=norm, cmap='viridis')
# plt.colorbar(im, ax=ax_cutout)  # Add colorbar to the current axes
# ax_cutout.grid(color='white', ls='solid')
# ax_cutout.set_xlabel('RA (J2000)')
# ax_cutout.set_ylabel('DEC (J2000)')
# ax_cutout.set_title('test af wise w3 cutout 1')
# ax_cutout.contour(cutout_data, levels=levels, colors='red', alpha=0.5)
