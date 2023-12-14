#%%
import pandas as pd
import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

################
### Get Flux ###
################

# load Scuba data
hdu_scuba = fits.open('C:/Users/olive/Dropbox/Group-1/Source Catalogs/maps/3C239_850_mf_cal_crop_snr.fits')
scuba850_data_cube = hdu_scuba[0].data  # This is a 3D cube
hdr_scuba = hdu_scuba[0].header
wcs_scuba = WCS(hdr_scuba)
df_scuba = pd.DataFrame(scuba850_data_cube[0,:,:])  # Pandas dataframe

# RAGERS coordinates txt file
positions_file = 'C:/Users/olive/Dropbox/5. semester/Special kursus/RAGERStxt_uden_obj36.txt'
txt_coord = pd.read_csv(positions_file, sep=" ")  # Pandas dataframe
obj_id = txt_coord['obj_id']

# WISE data
WISE = pd.read_csv('C:/Users/olive/Dropbox/5. semester/Special kursus/table_irsa_catalog_search_results.csv')
WISE_coord = WISE[['ra','dec']]

#####################
# Find target sources around RAGERS

# Input: WISE dataframe - RAGER ra and dec in deg - Other source ra and dec in deg - radius in arcsec
def ang_dist(df, RAGER_ra, RAGER_dec, source_ra, source_dec, radius):
    df['angular_distance'] = np.arccos(np.sin(RAGER_dec*np.pi/180)*np.sin(source_dec*np.pi/180)+np.cos(RAGER_dec*np.pi/180)*np.cos(source_dec*np.pi/180)*np.cos(RAGER_ra*np.pi/180 - source_ra*np.pi/180))*180/np.pi
    df['mask'] = df['angular_distance'] <= radius
    return df['mask'], df['angular_distance']

# Define coordinates
RAGER_ra = txt_coord['ra']  # RAGERS x_coord
RAGER_dec = txt_coord['dec']  # RAGERS y_coord
WISE_coord = WISE[['ra','dec']]

# Create n dataframes: Will contain data on sources within 15 arcsec of RAGERS
RAGERS = [['RAGER'+ str(_dummy)] for _dummy in range(1, 25)]

# Adds columns with names mask1 - mask24 including boolean if source is within target radius
# Adds columns with names ang_dist1 - ang_dist24 including boolean if source is within target radius
for i in range(0, len(txt_coord['ra'])):
    mask_id = ('mask' + str(i+1))
    ang_dist_id = ('ang_dist' + str(i+1))
    WISE[str(mask_id)] = ang_dist(WISE, RAGER_ra[i], RAGER_dec[i], WISE_coord['ra'], WISE_coord['dec'], 15.0/3600.0)[0]
    WISE[str(ang_dist_id)] = ang_dist(WISE, RAGER_ra[i], RAGER_dec[i], WISE_coord['ra'], WISE_coord['dec'], 15.0/3600.0)[1]

    # Loop to create a list of dataframes with sources within target radius 15 arcsec of RAGERS
    RAGERS[i] = pd.DataFrame(WISE[WISE.columns[0:36]][WISE[str(mask_id)]])
    RAGERS[i]['ang_dist'] = WISE[str(ang_dist_id)][WISE[str(mask_id)]]

# Count sources around each RAGER
n_source = np.array(np.zeros(len(txt_coord['ra'])))

for i in range(0, len(txt_coord['ra'])):
    if any(RAGERS[i]['w1mpro']):
        n_source[i] = np.array([len(RAGERS[i]['w1mpro'])])
    else:
        n_source[i] = 0

# Get source name for sources within 15 arcseconds
for i in range(0, len(txt_coord['ra'])):
    if any(RAGERS[i]['w1mpro']):
        WISE_names = RAGERS[i]['designation'].values[:]
        

# Print angular distances for each object
for i in range(len(txt_coord['ra'])):
    print(f"Object {obj_id[i]}:")
    print(RAGERS[i]['ang_dist'].values)
    print("\n")

# Create a list to store the flux values for each channel
flux_values = []

# Loop to calculate flux for WISE sources around RAGERS
for i in range(0, len(txt_coord['ra'])):
    m_vega_W1 = RAGERS[i]['w1mpro']
    m_vega_W2 = RAGERS[i]['w2mpro']
    m_vega_W3 = RAGERS[i]['w3mpro']
    m_vega_W4 = RAGERS[i]['w4mpro']

    # Constants
    Fv0_W1 = 309.540
    Fv0_W2 = 171.757
    Fv0_W3 = 31.678
    Fv0_W4 = 8.363

    # Check if there are any values in the arrays
    if not m_vega_W1.empty and not m_vega_W2.empty and not m_vega_W3.empty and not m_vega_W4.empty:
        # Equation 1 - Vega Magnitudes to Flux Density
        flux_W1 =  Fv0_W1 * 10**(-m_vega_W1 / 2.5)
        flux_W2 =  Fv0_W2 * 10**(-m_vega_W2 / 2.5)
        flux_W3 =  Fv0_W3 * 10**(-m_vega_W3 / 2.5)
        flux_W4 =  Fv0_W4 * 10**(-m_vega_W4 / 2.5)

        # Append the flux values to the list
        flux_values.append([obj_id[i], flux_W1.values[0], flux_W2.values[0], flux_W3.values[0], flux_W4.values[0]])

# Save the flux values to a text file
output_file = 'flux_values.txt'
with open(output_file, 'w') as f:
    f.write('obj_id\tflux_W1[r\"$\mu Jy$"]\tflux_W2[r\"$\mu Jy$"]\tflux_W3[r\"$\mu Jy$"]\tflux_W4[r\"$\mu Jy$"]\n')
    for entry in flux_values:
        f.write('\t'.join(map(str, entry)) + '\n')

print(f'Flux values saved to {output_file}')


# Save angular distances to a text file
output_ang_dist_file = 'angular_distances.txt'
with open(output_ang_dist_file, 'w') as f:
    f.write('obj_id\tangular_distance\n')
    for i in range(len(txt_coord['ra'])):
        obj_id_i = obj_id[i]
        ang_dist_values = RAGERS[i]['ang_dist'].values
        f.write(f"{obj_id_i}\t{', '.join(map(str, ang_dist_values))}\n")

print(f'Angular distances saved to {output_ang_dist_file}')
# %%
