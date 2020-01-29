import os
import sys
import numpy as np

import matplotlib as mpl
import matplotlib.cm as cm
import h5py

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

import flares_1.analysis.plt as fplt



DataFolder = f'/Users/stephenwilkins/research/simulations/FLARES/data/'

data = h5py.File(f'{DataFolder}/flares_noparticles.hdf5', 'r')


snaps = {'010_z005p000': 5.0, '009_z006p000': 6.0, '008_z007p000': 7.0, '007_z008p000': 8.0, '006_z009p000': 9.0, '005_z010p000': 10.0}

snap = '010_z005p000'

xlims = [-0.5,2.0]
ylims = [-0.5,2.0]



# --- all simulations (un-weighted)

all = False
if all:
    fig, ax = fplt.single()

    M = {}
    for f in ['U','V','J']:
        M[f] = np.array([])
        for sim in data.keys():
            M[f] = np.hstack((M[f], data[f'{sim}/{snap}/Galaxy/BPASS/ChabrierIMF/Luminosity/DustModelI/{f}']))

    UV = -2.5*np.log10(M['U']/M['V'])
    VJ = -2.5*np.log10(M['V']/M['J'])

    # ax.scatter(np.log10(mstar), np.log10(sfr), s=1, alpha = 0.25, c='k', zorder = 1)
    ax.hexbin(UV, VJ, gridsize = (25,25), bins = 'log', cmap='Greys', linewidths=0., mincnt = 1, extent = [*xlims, *ylims], alpha = 0.5, zorder = 2)

    ax.text(0.8, 0.9, f'$z={snaps[snap]}$', transform=ax.transAxes)

    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.set_ylabel(r'$U-V$')
    ax.set_xlabel(r'$V-J$')
    ax.grid(True)
    fig.savefig('figs/UVJ-all.pdf')



# --- color coded by sSFR

sSFR_color = True
if sSFR_color:
    fig, ax, cax = fplt.single_wcbar()

    # --- extract luminosities of all galaxies

    M = {}
    for f in ['U','V','J']:
        M[f] = np.array([])
        for sim in data.keys():
            M[f] = np.hstack((M[f], data[f'{sim}/{snap}/Galaxy/BPASS/ChabrierIMF/Luminosity/DustModelI/{f}']))

    # --- extract masses and sSFRs of all galaxies

    ssfr = np.array([])
    mstar = np.array([])
    for sim in data.keys():
        v = data[f'{sim}/{snap}/Galaxy']['SFR']['SFR_10'][:]/data[f'{sim}/{snap}/Galaxy']['Mstar_30'][:]
        ssfr = np.hstack((ssfr, v))
        mstar = np.hstack((mstar, data[f'{sim}/{snap}/Galaxy']['Mstar_30']))


    # --- calculate U-V and V-J color of galaxies

    UV = -2.5*np.log10(M['U']/M['V'])
    VJ = -2.5*np.log10(M['V']/M['J'])


    ssfr = np.log10(ssfr)-1 # --- correct SFR for duration (WILL NOT BE NEEDED IN FUTURE)
    ssfr += 9 # convert to Gyr^-1

    # --- apply mass cut

    s = np.log10(mstar)>9.

    # --- add scatter plot color-coded by sSFR

    norm = mpl.colors.Normalize(vmin=-0.5, vmax=1.5) # +1 gets rid of the yellow colour

    ax.scatter(VJ[s], UV[s], c = cm.coolwarm_r(norm(ssfr[s])), s=1, alpha = 0.5, zorder = 1)

    # --- UVJ selection region

    UVJ_region = [[-1., 0.75, 1.4, 1.4], [1.2, 1.2, 1.75, 2.0]]

    ax.fill_between(*UVJ_region,[xlims[1]]*len(UVJ_region[0]), color='k', alpha=0.1)
    ax.plot(*UVJ_region, c='k', lw=1)


    # --- add redshift label

    ax.text(0.8, 0.9, f'$z={snaps[snap]}$', transform=ax.transAxes)

    # --- add colorbar

    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cm.coolwarm_r), cax=cax)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'$\log_{10}({\rm sSFR}/{\rm Gyr^{-1})}$', fontsize=8)

    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.set_ylabel(r'$\rm U-V$')
    ax.set_xlabel(r'$\rm V-J$')
    ax.grid(True)


    fig.savefig('figs/UVJ-ssfr.pdf')


sSFR_color_intrinsic = True
if sSFR_color_intrinsic:
    fig, ax, cax = fplt.single_wcbar()

    # --- extract luminosities of all galaxies

    M = {}
    for f in ['U','V','J']:
        M[f] = np.array([])
        for sim in data.keys():
            M[f] = np.hstack((M[f], data[f'{sim}/{snap}/Galaxy/BPASS/ChabrierIMF/Luminosity/Intrinsic/{f}']))

    # --- extract masses and sSFRs of all galaxies

    ssfr = np.array([])
    mstar = np.array([])
    for sim in data.keys():
        v = data[f'{sim}/{snap}/Galaxy']['SFR']['SFR_10'][:]/data[f'{sim}/{snap}/Galaxy']['Mstar_30'][:]
        ssfr = np.hstack((ssfr, v))
        mstar = np.hstack((mstar, data[f'{sim}/{snap}/Galaxy']['Mstar_30']))


    # --- calculate U-V and V-J color of galaxies

    UV = -2.5*np.log10(M['U']/M['V'])
    VJ = -2.5*np.log10(M['V']/M['J'])


    ssfr = np.log10(ssfr)-1 # --- correct SFR for duration (WILL NOT BE NEEDED IN FUTURE)
    ssfr += 9 # convert to Gyr^-1

    # --- apply mass cut

    s = np.log10(mstar)>9.

    # --- add scatter plot color-coded by sSFR

    norm = mpl.colors.Normalize(vmin=-0.5, vmax=1.5) # +1 gets rid of the yellow colour

    ax.scatter(VJ[s], UV[s], c = cm.coolwarm_r(norm(ssfr[s])), s=1, alpha = 0.5, zorder = 1)

    # --- UVJ selection region

    UVJ_region = [[-1., 0.75, 1.4, 1.4], [1.2, 1.2, 1.75, 2.0]]

    ax.fill_between(*UVJ_region,[xlims[1]]*len(UVJ_region[0]), color='k', alpha=0.1)
    ax.plot(*UVJ_region, c='k', lw=1)


    # --- add redshift label

    ax.text(0.8, 0.9, f'$z={snaps[snap]}$', transform=ax.transAxes)

    # --- add colorbar

    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cm.coolwarm_r), cax=cax)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'$\log_{10}({\rm sSFR}/{\rm Gyr^{-1})}$', fontsize=8)

    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.set_ylabel(r'$\rm U-V$')
    ax.set_xlabel(r'$\rm V-J$')
    ax.grid(True)


    fig.savefig('figs/UVJ-ssfr-intrinsic.pdf')
