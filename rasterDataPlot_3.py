# -*- coding: utf-8 -*-


### import packages needed
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import gaussian_kde
import rasterstats
import rasterio as rio

import mpl_scatter_density
from matplotlib.colors import LinearSegmentedColormap

from numpy import polyfit

# Marius HUBER, 10/05/2023
# _____________________________________________________________________________________________________________________

### data import

aspect_fileName = 'aspect2.tif'
slope_fileName = 'slope2.tif'

with rio.open(aspect_fileName, 'r') as dset:
    aspect = dset.read()

with rio.open(slope_fileName, 'r') as dset:
    slope = dset.read()


# _____________________________________________________________________________________________________________________

### some initial data treatment

## transform to 1D-array

aspect = np.ravel(aspect)
slope = np.ravel(slope)

aspect = np.deg2rad(aspect)


## eliminate NaN-values

NaNvalue = [-9999.0]
mask1 = np.isin(aspect, NaNvalue)
mask2 = np.isin(slope, NaNvalue)

aspect = np.delete(aspect, mask2)
slope = np.delete(slope, mask2)

# print(aspect)
# print(slope)


## create second aspect array for alternative interpolation
# | create a copy of the aspect array that is shifted from the origin [°] to overcome the non-periodic
# | boundary condition

aspect_shift = aspect
slope_shift = slope

aspect_shift = aspect_shift + np.pi             # shift by π (half a circle)
# | bring values back to (0-2π)-range
aspect_shift = np.where(aspect_shift > 2*np.pi, aspect_shift - 2*np.pi, aspect_shift)


# _____________________________________________________________________________________________________________________

### interpolation step for density mapping
# | URL of web pages I used for help:
# | https://stackoverflow.com/questions/70646328/how-to-create-a-polar-density-plot-having-a-list-of-angles-values- \
# |      in-degrees-in
# | https://stackoverflow.com/questions/31051882/matplotlib-density-plot-in-polar-coordinates/31142319#31142319


## set the number of interpolation steps: choose only (multiples of 4) + 1
# | this needs to be done to avoid overlapping issues of the two independent interpolations
# | this can be improved in future versions

N = 101
# N = 41
# N = 13


## this is a factor that can be set for the gaussian_kde function
# | if the factor is higher the smoothing is less

interpFactor_divisor = 1
# interpFactor_divisor = 3


## compute the interpolation (2 times)

interp1 = gaussian_kde(np.vstack((aspect, slope)), bw_method='scott')
interp1.set_bandwidth(bw_method=interp1.factor / interpFactor_divisor)

interp2 = gaussian_kde(np.vstack((aspect_shift, slope_shift)), bw_method='scott')
interp2.set_bandwidth(bw_method=interp2.factor / interpFactor_divisor)


## define ranges for the density map, in which the interpolation should be displayed

# | aspect range for 360° in radians
aspect_ = np.linspace(0, 2*np.pi, N)
aspect_shift_ = np.linspace(0, 2*np.pi, N)

# | slope range, set upper bound in [°]
upperRbound = 60
slope_ = np.linspace(0, upperRbound, N)

## get the meshgrid
# | quote from source: @to compute colors, first get the meshgrid of angles and radii with shape N x N x 2"

mesh1 = np.stack(np.meshgrid(aspect_, slope_), 0)
mesh2 = np.stack(np.meshgrid(aspect_shift_, slope_), 0)

## reshaping of mesh (necessary?)
# | quote from source: "then, compute N x N matrix of colors using `interp`, \
# |     I think we need to reshape to accomodate Gaussian KDE API, but maybe this can be avoided"

Z1 = interp1(mesh1.reshape(2, -1)).reshape(N, N)
Z2 = interp2(mesh2.reshape(2, -1)).reshape(N, N)


## realign second aspect array to their original orientations

aspect_shift_ = aspect_shift_ - np.pi
# aspect_shift_ = np.where(aspect_shift_ < 0, aspect_shift_ + 2*np.pi, aspect_shift_)
print(aspect_shift_)

# _____________________________________________________________________________________________________________________

### perform a clipping on two independent density maps to visualize them on the same plot

## clipping mask for cutting off half a circle (180°)

mask3 = np.where(aspect_ < (1/2)*np.pi)
mask4 = np.where(aspect_ >= (3/2)*np.pi)
mask5 = np.hstack((mask3, mask4))


## performing the clipping on interpolation results (z-value arrays) and rearrange to a single new array

Z1_cut = np.delete(Z1, mask5, 1)
Z2_cut = np.delete(Z2, mask5, 1)

# print(Z1_cut.shape)
# print(Z2_cut.shape)
# print(Z1_cut)
# print(Z2_cut)

## further cutting the second interpolation results array to allow plotting on 0 - 2π range

Z2_cutCut1 = Z2_cut[:, 0:(N//4)]
Z2_cutCut2 = Z2_cut[:, (N//4):]

# print(Z2_cutCut1)
# print(Z2_cutCut2)
# print(Z2_cutCut1.shape)
# print(Z2_cutCut2.shape)


## Merge the interpolation results in a single new array

Z_combined = np.concatenate((Z2_cutCut2, Z1_cut, Z2_cutCut1, Z2_cutCut2[:, [0]]), axis=1)

# print(Z_combined.shape)
# print(Z_combined)


## Extract the maximum interpolation value for each column of the interpolation results array (represent maximum of
# | slope-angles)


maxLine_arg_y = Z_combined.argmax(axis=0)
# maxLine_y =

print(maxLine_arg_y)

maxLine_y = np.take(slope_, maxLine_arg_y)

print(maxLine_y)


# _____________________________________________________________________________________________________________________

### plotting

## here to set some plotting params

plt.rcParams.update({'figure.autolayout': True, 'axes.labelsize': 12, 'legend.fontsize': 14, 'xtick.labelsize': 12,
                     'ytick.labelsize': 12})
# plt.rcParams.update({'legend.numpoints': 1, 'font.size': 16, 'axes.labelsize': 16, 'xtick.major.pad': 10,
#                      'ytick.major.pad': 10, 'legend.fontsize': 16, 'xtick.labelsize': 12, 'ytick.labelsize': 12})
# lw=2
# ms=10
# the levels have to be set accordingly by trying (not so easy)
levels = np.linspace(0, 7e-3, 8)
alphaValue = 1


## figure 1: interpolation plot

fig1Name = 'aspectSlopeA_densityPlot'
fig1 = plt.figure(figsize=(20, 20))
ax1 = fig1.add_subplot(1, 1, 1)

ctf1 = ax1.contourf(aspect_, slope_, Z_combined, alpha=alphaValue, levels=levels, extend='max')   #, cmap='Greys')
plt.colorbar(ctf1)

ax1.plot(aspect_, maxLine_y)

ax1.set_ylim(0, 60)
ax1.set_xticks(np.deg2rad([0, 45, 90, 135, 180, 225, 270, 315, 360]))


### figure 2: try the mpl-density scatter plot

plt.style.use('dark_background')        # introduce a dark background style

fig2Name = 'densityScatter'
fig2 = plt.figure(figsize=(8, 40))
ax3 = fig2.add_subplot(1, 1, 1, projection='scatter_density')


## mpl-scatter-density plot

# "Viridis-like" colormap with black background

black_viridis = LinearSegmentedColormap.from_list('white_viridis', [
    (0, '#000000'),
    (1e-20, '#440053'),
    (0.2, '#404388'),
    (0.4, '#2a788e'),
    (0.6, '#21a784'),
    (0.8, '#78d151'),
    (1, '#fde624'),
], N=256)

# plotting

density = ax3.scatter_density(aspect, slope, vmin=0., vmax=6., cmap=black_viridis)
# fig2.colorbar(density, label='Number of points per pixel')
ax3.set_xlim(0, 2*np.pi)
ax3.set_ylim(0, 70)
ax3.set_xticks(np.deg2rad([0, 45, 90, 135, 180, 225, 270, 315, 360]))


# adding secondary ticks and labels and to axes, and a grid

def deg2rad(x):
    return x * np.pi / 180
def rad2deg(x):
    return x * 180 / np.pi

secax = ax3.secondary_xaxis('top', functions=(rad2deg, deg2rad))
secax.set_xticks([0, 45, 90, 135, 180, 225, 270, 315, 360])

secay = ax3.secondary_yaxis('right')
ax3.grid(axis='y')


# adding a polynomial regression line for maximum values of slope for interpolation

polyreg_maxLine = polyfit(aspect_, maxLine_y, deg=9)
polyreg_maxLine = np.poly1d(polyreg_maxLine)
polyreg_x=np.linspace(0, 2*np.pi, num=N)

ax3.plot(polyreg_x, polyreg_maxLine(polyreg_x), color='red', ls=':', linewidth=4)
ax3.scatter(aspect_, maxLine_y, color='red')


## export plots

# # figure 1
# fig1.savefig(fig1Name + '_N' + str(N) + '.eps', dpi=150, format='eps', transparent=False)
# fig1.savefig(fig1Name + '_N' + str(N) + '.png', dpi=150, format='png', transparent=False)
#
# # figure 2
# fig2.savefig(fig2Name + '_N' + str(N) + '.eps', dpi=300, format='eps', transparent=False)
# fig2.savefig(fig2Name + '_N' + str(N) + '.png', dpi=300, format='png', transparent=False)


# ----------------------------------------------------------------------------------------------
#### show
plt.show()
