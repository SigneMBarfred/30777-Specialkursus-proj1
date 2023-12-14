#%%
import pandas as pd
import numpy as np

# ... (your other imports)

# Load IRAC data
IRAC = pd.read_csv('C:/Users/olive/Dropbox/5. semester/Special kursus/IRAC_source_cat_3c239_400arcsec_radius.csv')

Flux_CH1_JY = IRAC['i1_f_ap1'] * 10**(-6)
Flux_CH2_JY = IRAC['i2_f_ap1'] * 10**(-6)
Flux_CH3_JY = IRAC['i3_f_ap1'] * 10**(-6)
Flux_CH4_JY = IRAC['i4_f_ap1'] * 10**(-6)

# Zero magnitude flux density
F0_3_6 = 280.9  # Jy
F0_4_5 = 179.7  # Jy
F0_5_8 = 115.0  # Jy
F0_8_0 = 64.9   # Jy

# Check for non-positive flux values
valid_flux_CH1 = Flux_CH1_JY > 0
valid_flux_CH2 = Flux_CH2_JY > 0
valid_flux_CH3 = Flux_CH3_JY > 0
valid_flux_CH4 = Flux_CH4_JY > 0

# Calculate magnitudes for valid flux values
m_3_6 = np.full_like(Flux_CH1_JY, np.nan)
m_4_5 = np.full_like(Flux_CH2_JY, np.nan)
m_5_8 = np.full_like(Flux_CH3_JY, np.nan)
m_8_0 = np.full_like(Flux_CH4_JY, np.nan)

m_3_6[valid_flux_CH1] = -0.1085736205e1 * np.log(Flux_CH1_JY[valid_flux_CH1] / F0_3_6)
m_4_5[valid_flux_CH2] = -0.1085736205e1 * np.log(Flux_CH2_JY[valid_flux_CH2] / F0_4_5)
m_5_8[valid_flux_CH3] = -0.1085736205e1 * np.log(Flux_CH3_JY[valid_flux_CH3] / F0_5_8)
m_8_0[valid_flux_CH4] = -0.1085736205e1 * np.log(Flux_CH4_JY[valid_flux_CH4] / F0_8_0)

# find highest magnitude in each band
mag = np.nanmax([m_3_6, m_4_5, m_5_8, m_8_0])

print(np.max(m_3_6[valid_flux_CH1]))
print(np.max(m_4_5[valid_flux_CH2]))
print(np.max(m_5_8[valid_flux_CH3]))
print(np.max(m_8_0[valid_flux_CH4]))


# %%
