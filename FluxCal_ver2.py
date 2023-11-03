
import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
import astropy.units as u
import math

# Load SCUBA data
hdu_scuba = fits.open('C:/Users/olive/Dropbox/Group-1/Source Catalogs/maps/3C239_850_mf_cal_crop_snr.fits')
scuba850_data_cube = hdu_scuba[0].data
hdr_scuba = hdu_scuba[0].header
wcs_scuba = WCS(hdr_scuba)

# Load RAGERS coordinates
positions_file = 'C:/Users/olive/Dropbox/5. semester/Special kursus/RAGERStxt_uden_obj36.txt'
txt_coord = pd.read_csv(positions_file, sep=" ")
obj_id = txt_coord['obj_id']

# Load WISE data
WISE = pd.read_csv('C:/Users/olive/Dropbox/5. semester/Special kursus/table_irsa_catalog_search_results.csv')

# Loop through RAGERS coordinates
for i in range(len(txt_coord)):
    x_coord = txt_coord['ra'][i]
    y_coord = txt_coord['dec'][i]
    target_radius = 15 / 3600  # 15 arcseconds
    
    # Calculate distance to WISE sources
    WISE['distance_to_RAGERS'] = np.sqrt((WISE['ra'] - x_coord) ** 2 + (WISE['dec'] - y_coord) ** 2)
    
    # Create a mask for sources within the target radius
    WISE[f'mask{i+1}'] = WISE['distance_to_RAGERS'] <= target_radius

    # Calculate flux for WISE sources within the target radius
    for band in ['W1', 'W2', 'W3', 'W4']:
        m_vega = WISE[f'{band.lower()}mpro']
        m_vega_error = WISE[f'{band.lower()}sigmpro']
        Fv0 = {'W1': 309.540, 'W2': 171.757, 'W3': 31.678, 'W4': 8.363}[band]
        
        # Calculate flux
        WISE[f'flux_{band}_RAGERS{i+1}'] = Fv0 * 10 ** (-m_vega / 2.5)
        
        # Calculate flux error
        WISE[f'flux_error_{band}_RAGERS{i+1}'] = Fv0 * 10 ** (-m_vega_error / 2.5)

# Count sources around each RAGERS
n_sources = [WISE[f'mask{i+1}'].sum() for i in range(len(txt_coord))]

# Access the calculated fluxes and errors for sources around each RAGERS
for i in range(len(txt_coord)):
    for band in ['W1', 'W2', 'W3', 'W4']:
        flux = WISE[f'flux_{band}_RAGERS{i+1}']
        flux_error = WISE[f'flux_error_{band}_RAGERS{i+1}']
        print(f'Flux and Error for RAGERS {i+1} and band {band}:')
        print(flux)
        print(flux_error)
