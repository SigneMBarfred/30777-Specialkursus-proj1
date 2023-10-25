#%%

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
fits_file_path = 'C:/Users/olive/Dropbox/5. semester/Special kursus/r10789376/ch1/pbcd/SPITZER_I1_10789376_0000_7_E8309859_maic.fits'
hdu_spitzer = fits.open(fits_file_path)
spitzer_data = hdu_spitzer[0].data
spitzer_hdr = hdu_spitzer[0].header
wcs_spitzer = WCS(spitzer_hdr)


import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.colors as colors
import os

# Define the base directory
base_dir = 'C:/Users/olive/Dropbox/5. semester/Special kursus/r10789376'

# Define the subfolder names
subfolders = ["ch1", "ch2", "ch3", "ch4"]

# Loop through the subfolders
for subfolder in subfolders:
    folder_path = os.path.join(base_dir, subfolder, 'pbcd')
    
    # Check if the folder exists
    if os.path.exists(folder_path):
        # List files in the folder and filter for files ending with "maic.fits"
        maic_files = [f for f in os.listdir(folder_path) if f.endswith('maic.fits')]
        
        if maic_files:
            for maic_file in maic_files:
                fits_file_path = os.path.join(folder_path, maic_file)
                
                # Open the FITS file and extract data and header
                hdu_spitzer = fits.open(fits_file_path)
                spitzer_data = hdu_spitzer[0].data
                spitzer_hdr = hdu_spitzer[0].header
                wcs_spitzer = WCS(spitzer_hdr)

                # The rest of your script for processing the data and creating cutouts
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
                plt.title(f'Spitzer Image {subfolder}')
                plt.xlabel('RA (J2000)')
                plt.ylabel('DEC (J2000)')

                # Read positions from your file
                positions_file = 'C:/Users/olive/Dropbox/5. semester/Special kursus/RAGERStxt_uden_obj36.txt'
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

                    # Check if the entire cutout image is completely white
                    if not np.all(np.isnan(cutout.data)):
                        plt.figure()
                        ax_cutout = plt.subplot(projection=cutout.wcs)
                        im = ax_cutout.imshow(cutout.data, norm=norm, cmap='Greys')
                        plt.colorbar(im, ax=ax_cutout)
                        ax_cutout.grid(color='white', ls='solid')
                        ax_cutout.set_xlabel('RA (J2000)')
                        ax_cutout.set_ylabel('DEC (J2000)')
                    
                        # Modify the title to include the current subfolder
                        ax_cutout.set_title(f'IRAC {subfolder} map of Obj_id {int(idno[count])} \n at {position}\n')

                        #cutout_filename = f'cutout_spitzer_{subfolder}_id_{int(idno[count])}.fits'
                        #fits.writeto(cutout_filename, cutout.data, overwrite=True)

                    count += 1

        else:
            print(f"No maic.fits files found in {subfolder}")
    else:
        print(f"{subfolder} folder not found")





