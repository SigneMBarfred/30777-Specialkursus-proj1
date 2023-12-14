#%%
import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

# Load IRAC data
IRAC = pd.read_csv('C:/Users/olive/Dropbox/5. semester/Special kursus/IRAC_source_cat_3c239_400arcsec_radius.csv')

# Unit definition
IRAC['i1_f_ap1'] = IRAC['i1_f_ap1'] * 10**(-6)
IRAC['i2_f_ap1'] = IRAC['i2_f_ap1'] * 10**(-6)
IRAC['i3_f_ap1'] = IRAC['i3_f_ap1'] * 10**(-6)
IRAC['i4_f_ap1'] = IRAC['i4_f_ap1'] * 10**(-6)

# RAGERS coordinates txt file
positions_file = 'C:/Users/olive/Dropbox/5. semester/Special kursus/RAGERStxt_uden_obj36.txt'
txt_coord = pd.read_csv(positions_file, sep=" ")
obj_id = txt_coord['obj_id']

# Define coordinates
RAGER_ra = txt_coord['ra']
RAGER_dec = txt_coord['dec']
IRAC_coord = IRAC[['ra', 'dec']]

# Input: IRAC dataframe - RAGER ra and dec in deg - Other source ra and dec in deg - radius in arcsec
def ang_dist_IRAC(df, RAGER_ra, RAGER_dec, source_ra, source_dec, radius):
    df['angular_distance'] = np.arccos(
        np.sin(RAGER_dec * np.pi / 180) * np.sin(source_dec * np.pi / 180) +
        np.cos(RAGER_dec * np.pi / 180) * np.cos(source_dec * np.pi / 180) *
        np.cos(RAGER_ra * np.pi / 180 - source_ra * np.pi / 180)
    ) * 180 / np.pi
    df['mask'] = df['angular_distance'] <= radius
    return df['mask'], df['angular_distance']

# Create n dataframes: Will contain data on sources within 15 arcsec of RAGERS
RAGERS_IRAC = [pd.DataFrame() for _dummy in range(len(txt_coord['ra']))]

# Loop to create list of dataframes with sources within target radius 15 arcsec of RAGERS for IRAC
for i in range(len(txt_coord['ra'])):
    mask_id = f'mask{i + 1}'
    ang_dist_id = f'ang_dist{i + 1}'
    
    # Call the ang_dist_IRAC function for the current RAGER
    mask, ang_dist_values = ang_dist_IRAC(IRAC, RAGER_ra[i], RAGER_dec[i], IRAC['ra'], IRAC['dec'], 15.0 / 3600.0)
    
    # Apply the mask and store angular distances in the IRAC DataFrame
    IRAC[mask_id] = mask
    IRAC[ang_dist_id] = ang_dist_values
    
    # Add sources to RAGERS DataFrame
    RAGERS_IRAC[i] = IRAC.loc[mask, IRAC.columns[0:56]]
    


# Count sources around each RAGER for IRAC
n_source_IRAC = np.zeros(len(txt_coord['ra']), dtype=int)

for i in range(len(txt_coord['ra'])):
    n_source_IRAC[i] = len(RAGERS_IRAC[i])
    print(f'RAGER{i + 1}: obj_id {obj_id[i]} has {n_source_IRAC[i]} sources within 15 arcseconds')

# Get source name for sources within 15 arcseconds for IRAC
for i in range(len(txt_coord['ra'])):
    if not RAGERS_IRAC[i].empty:
        IRAC_names = RAGERS_IRAC[i]['objid'].values
        print(f'RAGER{i + 1}: obj_id {obj_id[i]} is near {IRAC_names} ')

# Loop to calculate flux for IRAC sources around RAGERS
flux_values_IRAC = []

for i in range(len(txt_coord['ra'])):
    flux_values_IRAC.append([obj_id[i]] + RAGERS_IRAC[i][['i1_f_ap1', 'i2_f_ap1', 'i3_f_ap1', 'i4_f_ap1']].values.flatten().tolist())


# Save the flux values to a text file
output_file_IRAC = 'flux_values_IRAC_test.txt'
with open(output_file_IRAC, 'w') as f:
    f.write('obj_id\tflux_I1[r\"$\mu Jy$"]\tflux_I2[r\"$\mu Jy$"]\tflux_I3[r\"$\mu Jy$"]\tflux_I4[r\"$\mu Jy$"]\n')
    for entry in flux_values_IRAC:
        f.write('\t'.join(map(str, entry)) + '\n')

print(f'Flux values saved to {output_file_IRAC}')

# Save angular distances to a text file
output_ang_dist_file = 'angular_distances_IRAC.txt'
with open(output_ang_dist_file, 'w') as f:
    f.write('obj_id\tangular_distance\n')
    for i in range(len(txt_coord['ra'])):
        obj_id_i = obj_id[i]
        if obj_id_i in obj_ids_of_interest and not RAGERS_IRAC[i].empty:
            mask, ang_dist_values = ang_dist_IRAC(IRAC, RAGER_ra[i], RAGER_dec[i], IRAC_coord['ra'], IRAC_coord['dec'], 15.0 / 3600.0)
            angular_distances = ang_dist_values[mask]
            print(f"Object {obj_id_i}:")
            print(f"Angular distances: {angular_distances}")
            f.write(f"{obj_id_i}\t{', '.join(map(str, angular_distances))}\n")

print(f'Angular distances saved to {output_ang_dist_file}')




