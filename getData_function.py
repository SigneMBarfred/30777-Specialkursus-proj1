# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 10:09:06 2023

@author: Nikolaj Lange Dons
"""

from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import numpy as np
import os

# Input
# Define path to folder with data files
# folder_path = 'C:/Users/Nikolaj Lange Dons/OneDrive - Danmarks Tekniske Universitet/Dokumenter/5 semester/Space special kursus/phot/ch1/bcd/'
# and
# Data format f.eks bcd (string)

# output
# matrix with data from file

def getData(folder_path, format):
    
    # Identify file formats in folder
    folder_bcd = os.listdir(folder_path)

    format_file = (format + '_file')

    # sorter dataformat ind til matrix
    format_file = []

    for filename in folder_bcd:
        if filename[32:-5] == format:
            bcd_load = get_pkg_data_filename(folder_path + filename)
            format_file.append(bcd_load)

    # Loading data
    data = format
    data = []
    
    for i in range(0, np.size(format_file, 0)):
        format_data = fits.getdata(format_file[i], ext=0)
        data.append(format_data)
        
    return data

# Example

# folder_path = 'C:/Users/Nikolaj Lange Dons/OneDrive - Danmarks Tekniske Universitet/Dokumenter/5 semester/Space special kursus/phot/ch1/bcd/'
# format = 'bcd'
# bcd = getData(folder_path, format)


