#!/usr/bin/python
#
# Lacerda@Orotour - 24/Nov/2015
#
import sys
import h5py
import numpy as np
import CALIFAUtils as C
import matplotlib as mpl
from matplotlib import pyplot as plt
from CALIFAUtils.scripts import get_h5_data_masked
from CALIFAUtils.plots import plotOLSbisectorAxis

mpl.rcParams['font.size'] = 12
mpl.rcParams['axes.labelsize'] = 10
mpl.rcParams['axes.titlesize'] = 12
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

h5file = sys.argv[1]
h5 = h5py.File(h5file, 'r')

SFR_int = get_h5_data_masked(h5, 'SFR_int__Tg', **dict(dtype = np.float_))
SFR_Ha_int = get_h5_data_masked(h5, 'SFR_Ha_int__Tg', **dict(dtype = np.float_))
SFR_Ha_int_masked = get_h5_data_masked(h5, 'SFR_Ha_int_masked__Tg', **dict(dtype = np.float_))
O3N2_int_masked = get_h5_data_masked(h5, 'O3N2M13_int_masked__Tg', **dict(dtype = np.float_))
O3N2_int = np.array(h5['data/O3N2M13_int__g'].value, dtype = np.float_)
O3N2_int = np.ma.masked_array(O3N2_int, mask = np.isnan(O3N2_int))

sc_kwargs = dict(marker = 'o', s = 50, edgecolor = 'none', label = '')
kwargs_ols_plot = dict(c = 'r', ls = '--', lw = 2, label = 'OLS')
kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 12, rms = True, text = True, kwargs_plot = kwargs_ols_plot)

f = plt.figure()
grid_shape = (1, 2)
ax1 = plt.subplot2grid(grid_shape, loc = (0, 0))
ax2 = plt.subplot2grid(grid_shape, loc = (0, 1))
xm1, ym1 = C.ma_mask_xyz(x = O3N2_int, y = np.ma.log10(SFR_int[1]) - np.ma.log10(SFR_Ha_int[1]))
ax1.scatter(xm1, ym1, **sc_kwargs)
_ = plotOLSbisectorAxis(ax1, xm1.compressed(), ym1.compressed(), **kwargs_ols)
ax1.set_xlabel(r'12 + $\log$ O/H M13')
ax1.set_ylabel(r'$\Delta$SFR (syn - $H\alpha$)')
sc = ax2.scatter(np.ma.log10(SFR_int[1]), np.ma.log10(SFR_Ha_int[1]), c = O3N2_int, cmap = 'viridis', **sc_kwargs)
ax2.set_xlabel(r'$\log$ SFR(syn) [$M_\odot\ yr^{-1}$]')
ax2.set_ylabel(r'$\log$ SFR(H$\alpha$) [$M_\odot\ yr^{-1}$]')
f.colorbar(sc)
f.savefig('Zneb_SFR.pdf')
f.subplots_adjust(bottom = 0.15, hspace = 0.15, wspace = 0., right = 0.95, left = 0.10)

f = plt.figure()
grid_shape = (1, 2)
ax1 = plt.subplot2grid(grid_shape, loc = (0, 0))
ax2 = plt.subplot2grid(grid_shape, loc = (0, 1))
xm2, ym2 = C.ma_mask_xyz(x = O3N2_int_masked[1], y = np.ma.log10(SFR_int[1]) - np.ma.log10(SFR_Ha_int_masked[1]))
ax1.scatter(xm2, ym2, **sc_kwargs)
# Just putting the same scale relation in plot 2
_ = plotOLSbisectorAxis(ax1, xm1.compressed(), ym1.compressed(), **kwargs_ols)
ax1.set_xlabel(r'12 + $\log$ O/H M13')
ax1.set_ylabel(r'$\Delta$SFR (syn - $H\alpha$)')
sc = ax2.scatter(np.ma.log10(SFR_int[1]), np.ma.log10(SFR_Ha_int_masked[1]), c = O3N2_int_masked[1], cmap = 'viridis', **sc_kwargs)
ax2.set_xlabel(r'$\log$ SFR(syn) [$M_\odot\ yr^{-1}$]')
ax2.set_ylabel(r'$\log$ SFR(H$\alpha$) [$M_\odot\ yr^{-1}$]')
f.colorbar(sc)
f.savefig('Zneb_SFR_masked.pdf')
