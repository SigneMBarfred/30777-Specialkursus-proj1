# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 14:39:13 2023

@author: signe
"""
import matplotlib.colors
from astropy.io import fits
from matplotlib import pyplot as plt, cm
from matplotlib import colors
from astropy.visualization import astropy_mpl_style
from astropy.utils.data import get_pkg_data_filename

image_file = get_pkg_data_filename('tutorials/FITS-images/HorseHead.fits')

plt.style.use(astropy_mpl_style)

#loading data
#hdu = fits.open('C:/Users/signe/Documents/DTU/Specialkursus_radiogal/r10760192/ch1/bcd/SPITZER_I1_10760192_0000_0000_6_bcd.fits')


image_file = get_pkg_data_filename('r10760192/ch1/bcd/SPITZER_I1_10760192_0000_0000_6_bcd.fits')
fits.info(image_file)
image_data = fits.getdata(image_file, ext=0)


colors = [[1,.5,.5], [1,0,0], [0,1,0], [0,0,1], [0,0,0], [0, .5, .5]]
cm = matplotlib.colors.ListedColormap(colors)
norm = matplotlib.colors.BoundaryNorm([0,1,2,4,6,8,12], cm.N)


plt.figure()
im = plt.imshow(image_data, cmap=cm, norm=norm)
#im.set_clim()
plt.colorbar(im)