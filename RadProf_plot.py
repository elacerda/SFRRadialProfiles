#!/usr/bin/python
#
# Lacerda@Corrego - 21/Nov/2015
#
import sys
import h5py
import numpy as np
import matplotlib as mpl
import CALIFAUtils as C
from matplotlib import pyplot as plt
from CALIFAUtils.plots import plot_text_ax
from matplotlib.pyplot import MultipleLocator
from CALIFAUtils.scripts import get_h5_data_masked

if __name__ == '__main__':
    try:
        h5file = sys.argv[1]
    except IndexError:
        print 'usage: %s HDF5FILE' % (sys.argv[0])
        exit(1)
    
    h5 = h5py.File(h5file, 'r')
    
    RbinCenter__r = h5['HEADER/RbinCenter__r'].value
    rbinstep = h5['HEADER/RbinStep'].value
    tSF__T = h5[ 'HEADER/tSF__T'].value 
    NRbins = h5['HEADER/N_rbins'].value
    N_gals = h5['HEADER/N_gals']
    N_T = h5['HEADER/N_T'].value
    gals = h5['HEADER/gals'].value
    morf__g = get_h5_data_masked(h5, 'morf__g', **dict(dtype = np.float_))
    aSFRSD__Trg = get_h5_data_masked(h5, 'aSFRSD__Trg', **dict(dtype = np.float_))
    aSFRSD_Ha__Trg = get_h5_data_masked(h5, 'aSFRSD_Ha__Trg', **dict(dtype = np.float_))
    N_zones_notmasked__Tg = get_h5_data_masked(h5, 'N_zones_notmasked__Tg', **dict(dtype = np.float_))

    colortipo = ['brown', 'red', 'orange', 'green', '#00D0C9', '#0076C9', 'blue']
    Ntype = len(colortipo)
    mtypes = [ 0, 1, 2, 3, 4, 5, 6 ]
    mtype_labels = [ 'E', 'S0', 'Sa', 'Sb', 'Sbc', 'Sc', 'Sd' ]
    halfbinstep = np.diff(mtypes)[0]/2.
    tickpos = np.linspace(mtypes[0] + halfbinstep, mtypes[-1] - halfbinstep, Ntype)
            
    aSFRSD__Ttr = np.empty((N_T, Ntype, NRbins), dtype = np.float_)
    aSFRSD_Ha__Ttr = np.empty((N_T, Ntype, NRbins), dtype = np.float_)
    nbinelem_aSFRSD__Ttr = np.empty((N_T, Ntype, NRbins), dtype = np.int_)
    nbinelem_aSFRSD_Ha__Ttr = np.empty((N_T, Ntype, NRbins), dtype = np.int_)
    perc_aSFRSD__pTtr = np.ma.masked_all((3, N_T, Ntype, NRbins), dtype = np.float_)
    perc_aSFRSD_Ha__pTtr = np.ma.masked_all((3, N_T, Ntype, NRbins), dtype = np.float_)
    
    Ngals__t = np.empty((Ntype), dtype = np.int_)
    
    for it in mtypes:
        mask_morf__g = (morf__g == it)
        Ngals__t[it] = np.sum(mask_morf__g)
        for iT in range(N_T):
    
            # G means masked g
            aSFRSD__rG = aSFRSD__Trg[iT][:, mask_morf__g]
            aSFRSD_Ha__rG = aSFRSD_Ha__Trg[iT][:, mask_morf__g]
    
            for ir in np.where(RbinCenter__r <= 2 + rbinstep)[0]:
                for ig in xrange(Ngals__t[it]):
                    if aSFRSD__rG[ir].mask[ig]:
                        aSFRSD__rG[ir, ig] = 0
                    if aSFRSD_Ha__rG[ir].mask[ig]:    
                        aSFRSD_Ha__rG[ir, ig] = 0.
                    
            aSFRSD__Ttr[iT][it] = aSFRSD__rG.mean(axis = 1)
            nbinelem_aSFRSD__Ttr[iT][it] = aSFRSD__rG.count(axis = 1)
            aSFRSD_Ha__Ttr[iT][it] = aSFRSD_Ha__rG.mean(axis = 1)
            nbinelem_aSFRSD_Ha__Ttr[iT][it] = aSFRSD_Ha__rG.count(axis = 1)
            
            for ir in xrange(NRbins):
                if nbinelem_aSFRSD__Ttr[iT, it, ir] > 2:
                    perc_aSFRSD__pTtr[:, iT, it, ir] = np.percentile(aSFRSD__rG[ir].compressed(), [16, 50, 84])
    
                if nbinelem_aSFRSD_Ha__Ttr[iT, it, ir] > 2:
                    perc_aSFRSD_Ha__pTtr[:, iT, it, ir] = np.percentile(aSFRSD__rG[ir].compressed(), [16, 50, 84])

    for iT in xrange(N_T):    
        ylim = (-5, 0)
        f = plt.figure()
        NRows = 2
        NCols = 8
        page_size_inches = [12, 8]
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
        cmap = mpl.colors.ListedColormap(colortipo)
        suptitle = 'Ngals:%d  Nzones:%d  NRbins:%d  tSF:%.2f Myr' % ( np.sum(morf__g >= 0), N_zones_notmasked__Tg[iT].sum(), aSFRSD__Trg[iT].count(), tSF__T[iT]/1e6)
        try:
            suptitle = ('%s  %' % (suptitle, sys.argv[3]))
        f.suptitle(suptitle, fontsize = 15) 
        
        ax = plt.subplot2grid(grid_shape, loc = (0, 0))
        for it in mtypes:
            ax.plot(RbinCenter__r, np.ma.log10(aSFRSD__Ttr[iT][it] * 1e6), c = colortipo[it], label = mtype_labels[it])
        ax.set_ylabel(r'$\log\ \Sigma_{SFR}^\star(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$')
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_minor_locator(MultipleLocator(.2))
        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_minor_locator(MultipleLocator(.2))
        plt.setp(ax.get_xticklabels(), visible = False)
        ax.grid(which = 'major')
        ax.set_xlim(0,3)
        ax.set_ylim(ylim)
        
        ax = plt.subplot2grid(grid_shape, loc = (1, 0))
        for it in mtypes:
            ax.plot(RbinCenter__r, np.ma.log10(aSFRSD_Ha__Ttr[iT][it] * 1e6), c = colortipo[it], label = mtype_labels[it])
        ax.set_ylabel(r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$')
        ax.set_xlabel(r'R [HLR]')
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_minor_locator(MultipleLocator(.2))
        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_minor_locator(MultipleLocator(.2))
        ax.grid(which = 'major')
        ax.set_xlim(0,3)
        ax.set_ylim(ylim)
        
        for it in mtypes:
            ax = plt.subplot2grid(grid_shape, loc = (0, it + 1))
            ax.plot(RbinCenter__r, np.ma.log10(aSFRSD__Ttr[iT][it] * 1e6), c = colortipo[it], label = mtype_labels[it])
            ax.plot(RbinCenter__r, np.ma.log10(perc_aSFRSD__pTtr[0, iT, it, :] * 1e6), 'k--')
            ax.plot(RbinCenter__r, np.ma.log10(perc_aSFRSD__pTtr[1, iT, it, :] * 1e6), 'k--')
            ax.plot(RbinCenter__r, np.ma.log10(perc_aSFRSD__pTtr[2, iT, it, :] * 1e6), 'k--')
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(.2))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(.2))
            plt.setp(ax.get_xticklabels(), visible = False)
            plt.setp(ax.get_yticklabels(), visible = False)
            plot_text_ax(ax, '%s (%d)' % (mtype_labels[it], Ngals__t[it]), 0.95, 0.95, 15, 'top', 'right', colortipo[it])
            ax.grid(which = 'major')
            ax.set_xlim(0,3)
            ax.set_ylim(ylim)
    
            ax = plt.subplot2grid(grid_shape, loc = (1, it + 1))
            ax.plot(RbinCenter__r, np.ma.log10(aSFRSD_Ha__Ttr[iT][it] * 1e6), c = colortipo[it], label = mtype_labels[it])
            ax.plot(RbinCenter__r, np.ma.log10(perc_aSFRSD_Ha__pTtr[0, iT, it, :] * 1e6), 'k--')
            ax.plot(RbinCenter__r, np.ma.log10(perc_aSFRSD_Ha__pTtr[1, iT, it, :] * 1e6), 'k--')
            ax.plot(RbinCenter__r, np.ma.log10(perc_aSFRSD_Ha__pTtr[2, iT, it, :] * 1e6), 'k--')
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(.2))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(.2))
            plt.setp(ax.get_xticklabels(), visible = False)
            plt.setp(ax.get_yticklabels(), visible = False)
            ax.grid(which = 'major')
            ax.set_xlim(0,3)
            ax.set_ylim(ylim)
            
            C.debug_var(
                True, 
                tipoM = mtypes[it], 
                Nelemperbin = nbinelem_aSFRSD__Ttr[iT][it],
                Nelemperbin_Ha = nbinelem_aSFRSD_Ha__Ttr[iT][it],
            )
        
        f.subplots_adjust(bottom = 0.15, hspace = 0, wspace = 0, right = 0.95, left = 0.07)

        try:
            output = '%s_%.2fMyrs.pdf' % (sys.argv[2], tSF__T[iT]/1e6)
        except IndexError:
            output = 'SFRSDRadialProfiles_%.2fMyrs.pdf' % tSF__T[iT]/1e6
            print 'no output name, using: %s' % (sys.argv[0])
            
        f.savefig(output)
