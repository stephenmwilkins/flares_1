import os
import sys
import numpy as np

import h5py

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

import flares_1.analysis.plt as fplt



DataFolder = f'/Users/stephenwilkins/research/simulations/FLARES/data/'

data = h5py.File(f'{DataFolder}/flares_noparticles.hdf5', 'r')


snaps = {'010_z005p000': 5.0, '009_z006p000': 6.0, '008_z007p000': 7.0, '007_z008p000': 8.0, '006_z009p000': 9.0, '005_z010p000': 10.0}

snap = '010_z005p000'

xlims = [8, 11.5]
ylims = [-1.5, 2.]


# --- coloured by simulation

by_sim = False
if by_sim:
    fig, ax = fplt.single()
    for sim in data.keys():

        sfr = data[f'{sim}/{snap}/Galaxy']['SFR']['SFR_10']
        mstar = data[f'{sim}/{snap}/Galaxy']['Mstar_30']

        ax.scatter(np.log10(mstar), np.log10(sfr)-np.log10(mstar), label = sim, s=5, alpha = 0.25)

    # ax.set_ylim([27, 32])
    # ax.set_xlim([7, 12])

    ax.set_ylabel(r'$\log_{10}({\rm SFR}/{\rm yr^{-1})}$')
    ax.set_xlabel(r'$\log_{10}(M_{\star}/{\rm M_{\odot}})$')
    ax.grid(True)
    ax.legend(fontsize=6)
    fig.savefig('figs/mass_sfr-sims.pdf')



# --- all simulations (un-weighted)

all = True
if all:
    fig, ax = fplt.single()

    sfr = np.array([])
    mstar = np.array([])

    for sim in data.keys():
        sfr = np.hstack((sfr, data[f'{sim}/{snap}/Galaxy']['SFR']['SFR_10']))
        mstar = np.hstack((mstar, data[f'{sim}/{snap}/Galaxy']['Mstar_30']))


    sfr[sfr<=0] = 1E-10
    sfr = np.log10(sfr)
    mstar = np.log10(mstar)

    ssfr = sfr-mstar-1+9

    # ax.scatter(np.log10(mstar), np.log10(sfr), s=1, alpha = 0.25, c='k', zorder = 1)
    ax.hexbin(mstar, ssfr, gridsize = (25,25), bins = 'log', cmap='Greys', linewidths=0., mincnt = 1, extent = [*xlims, *ylims], alpha = 0.5, zorder = 2)

    bins = np.arange(*xlims, 0.2)
    P16, P50, P84 = fplt.average_line(mstar, ssfr,bins)

    ax.fill_between(bins,P84,P16,color='k', alpha=0.15)
    ax.plot(bins,P50,ls='-',c='k', alpha=1.0, lw=1)

    ax.text(0.8, 0.9, f'$z={snaps[snap]}$', transform=ax.transAxes)

    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.set_ylabel(r'$\log_{10}({\rm sSFR}/{\rm Gyr^{-1})}$')
    ax.set_xlabel(r'$\log_{10}(M_{\star}/{\rm M_{\odot}})$')
    ax.grid(True)
    fig.savefig('figs/mass_sfr-all.pdf')



# --- all simulations (un-weighted)

redshift_evolution = True
if redshift_evolution:

    fig, ax = fplt.single()

    for snap, z in snaps.items():

        print('-'*5, snap)

        sfr = np.array([])
        mstar = np.array([])

        for sim in data.keys():
            sfr = np.hstack((sfr, data[f'{sim}/{snap}/Galaxy']['SFR']['SFR_10']))
            mstar = np.hstack((mstar, data[f'{sim}/{snap}/Galaxy']['Mstar_30']))

        sfr[sfr<=0] = 1E-10
        sfr = np.log10(sfr)
        mstar = np.log10(mstar)
        ssfr = sfr-mstar-1+9

        print(f'N: {sfr.shape[0]}')
        print(f'max SFR: {np.max(sfr):.2f}')
        print(f'max M*: {np.max(mstar):.2f}')

        bins = np.arange(*xlims, 0.2)
        P16, P50, P84 = fplt.average_line(mstar, ssfr,bins)

        ax.fill_between(bins,P84,P16, color=fplt.c_z(z), alpha=0.1)
        ax.plot(bins,P50,ls='-', color=fplt.c_z(z), alpha=1.0, lw=1, label = f'$z={z}$')



    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.set_ylabel(r'$\log_{10}({\rm sSFR}/{\rm Gyr^{-1})}$')
    ax.set_xlabel(r'$\log_{10}(M_{\star}/{\rm M_{\odot}})$')
    ax.legend()
    ax.grid(True)
    fig.savefig('figs/mass_sfr-z.pdf')
