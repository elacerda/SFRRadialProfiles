#!/usr/bin/python
#
# Lacerda@Granada - 14/Dec/2015 
#
import sys
import numpy as np
import matplotlib as mpl
from scipy import stats as st
from matplotlib import pyplot as plt
from CALIFAUtils.plots import plot_text_ax
from CALIFAUtils.scripts import OLS_bisector, ma_mask_xyz
from matplotlib.ticker import MultipleLocator
from CALIFAUtils.plots import plotOLSbisectorAxis

mpl.rcParams['font.size'] = 12
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['axes.titlesize'] = 14
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

if __name__ == '__main__':
    import h5py
    from CALIFAUtils.scripts import get_h5_data_masked
    baseCode = 'Bzca6e'
    #baseCode = 'Bgsd6e'
    
    h5 = h5py.File('/Users/lacerda/dev/astro/SFRRadialProfiles/SFRSDRadialProfiles_minpopx0.05_%s.h5' % baseCode)
    
    N_T = h5['HEADER/N_T'].value.astype(np.int)
    tSF__T = h5['HEADER/tSF__T'].value.astype(np.float_)
    SFR__Tg = [ get_h5_data_masked(h5, 'SFR__Tg/%d'%iT, **dict(dtype = np.float_)) for iT in xrange(N_T) ]
    SFRSD__Tg = [ get_h5_data_masked(h5, 'SFRSD__Tg/%d'%iT, **dict(dtype = np.float_)) for iT in xrange(N_T) ]
    SFR_Ha__Tg = [ get_h5_data_masked(h5, 'SFR_Ha__Tg/%d'%iT, **dict(dtype = np.float_)) for iT in xrange(N_T) ]
    SFRSD_Ha__Tg = [ get_h5_data_masked(h5, 'SFRSD_Ha__Tg/%d'%iT, **dict(dtype = np.float_)) for iT in xrange(N_T) ]
    
    f = plt.figure()
    page_size_inches = [10, 8]
    f.set_size_inches(page_size_inches)
    f.suptitle(r'%s  %d gals  minpopx 5%%' % (baseCode, h5['HEADER/N_gals'].value.astype(np.int)))
    default_sc_kwargs = dict(marker = 'o', s = 10, alpha = 0.8, edgecolor = 'none', label = '')
    default_ols_plot_kwargs = dict(c = 'r', ls = '--', lw = 2, label = 'OLS')
    default_ols_kwargs = dict(c = 'r', pos_x = 0.98, pos_y = 0.01, fs = 12, rms = True, text = True, kwargs_plot = default_ols_plot_kwargs)

    grid_shape = (2, 2)
    xlabel = r'$\log\ SFR^\star(t_\star)\ [M_\odot yr^{-1}]$'
    ylabel = r'$\log\ SFR^{neb}\ [M_\odot yr^{-1}]$'
    xSDlabel = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} pc^{-2}]$'
    ySDlabel = r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} pc^{-2}]$'
    
    #SFR
    iT = 1
    ax = plt.subplot2grid(grid_shape, loc = (0, 0))
    ax.set_title(r'$t_\star$ = %.2f Myrs' % (tSF__T[iT] / 1e6))
    xm, ym = ma_mask_xyz(np.ma.log10(SFR__Tg[iT]), np.ma.log10(SFR_Ha__Tg[iT]))
    sc = ax.scatter(xm, ym, **default_sc_kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(xlabel)
    ax.set_xlim([-6.5, 0])
    ax.set_ylim([-6.5, 0])
    a, b, sa, sb = plotOLSbisectorAxis(ax, xm, ym)
    ax.plot(ax.get_xlim(), ax.get_xlim(), 'k--')
    ax = plt.subplot2grid(grid_shape, loc = (0, 1))
    ax.set_title(r'$t_\star$ = %.2f Myrs' % (tSF__T[iT] / 1e6))
    xm, ym = ma_mask_xyz(np.ma.log10(SFRSD__Tg[iT] * 1e6), np.ma.log10(SFRSD_Ha__Tg[iT] * 1e6))
    sc = ax.scatter(xm, ym, **default_sc_kwargs)
    ax.set_xlabel(xSDlabel)
    ax.set_ylabel(xSDlabel)
    ax.set_xlim([-5, 1])
    ax.set_ylim([-5, 1])
    a, b, sa, sb = plotOLSbisectorAxis(ax, xm, ym)
    ax.plot(ax.get_xlim(), ax.get_xlim(), 'k--')

    iT = 2
    ax = plt.subplot2grid(grid_shape, loc = (1, 0))
    ax.set_title(r'$t_\star$ = %.2f Myrs' % (tSF__T[iT] / 1e6))
    xm, ym = ma_mask_xyz(np.ma.log10(SFR__Tg[iT]), np.ma.log10(SFR_Ha__Tg[iT]))
    sc = ax.scatter(xm, ym, **default_sc_kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(xlabel)
    ax.set_xlim([-6.5, 0])
    ax.set_ylim([-6.5, 0])
    a, b, sa, sb = plotOLSbisectorAxis(ax, xm, ym)
    ax.plot(ax.get_xlim(), ax.get_xlim(), 'k--')
    ax = plt.subplot2grid(grid_shape, loc = (1, 1))
    ax.set_title(r'$t_\star$ = %.2f Myrs' % (tSF__T[iT] / 1e6))
    xm, ym = ma_mask_xyz(np.ma.log10(SFRSD__Tg[iT] * 1e6), np.ma.log10(SFRSD_Ha__Tg[iT] * 1e6))
    sc = ax.scatter(xm, ym, **default_sc_kwargs)
    ax.set_xlabel(xSDlabel)
    ax.set_ylabel(xSDlabel)
    ax.set_xlim([-5, 1])
    ax.set_ylim([-5, 1])
    a, b, sa, sb = plotOLSbisectorAxis(ax, xm, ym)
    ax.plot(ax.get_xlim(), ax.get_xlim(), 'k--')
    
    f.subplots_adjust(hspace = 0.4)#, wspace = 0.2)
    f.savefig('SFR_%s.png' % baseCode)
    sys.exit(1)
    
    