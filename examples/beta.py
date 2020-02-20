import os
import sys
import numpy as np

import matplotlib as mpl
import matplotlib.cm as cm
import h5py

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

import flares_1.analysis.plt as fplt


Bouwens = {'log10LFUV' : np.array([29.33819871, 29.13819871, 28.93819871, 28.73819871, 28.53819871,
       28.33819871, 28.13819871, 27.93819871, 27.73819871, 27.43819871]), 'beta': np.array([-1.55, -1.58, -1.74, -1.9 , -1.9 , -2.22, -2.26, -2.19, -2.4 , -2.24])}



DataFolder = f'/Users/stephenwilkins/research/simulations/FLARES/data/'

data = h5py.File(f'{DataFolder}/flares_noparticles.hdf5', 'r')


snaps = {'010_z005p000': 5.0, '009_z006p000': 6.0, '008_z007p000': 7.0, '007_z008p000': 8.0, '006_z009p000': 9.0, '005_z010p000': 10.0}

snap = '010_z005p000'






# --- all simulations (un-weighted)

beta_luminosity = False
if beta_luminosity:

    xlims = [28, 30.5]
    ylims = [-2.75,-1.0]
    fig, ax = fplt.single()

    L = {}
    for f in ['FUV','NUV']:
        L[f] = np.array([])
        for sim in data.keys():
            L[f] = np.hstack((L[f], data[f'{sim}/{snap}/Galaxy/BPASS/ChabrierIMF/Luminosity/DustModelI/{f}']))


    beta = np.log10(L['FUV']/L['NUV'])/np.log10(1500./2500.) - 2.0

    # ax.scatter(np.log10(mstar), np.log10(sfr), s=1, alpha = 0.25, c='k', zorder = 1)
    ax.hexbin(np.log10(L['FUV']), beta, gridsize = (25,25), bins = 'log', cmap='Greys', linewidths=0., mincnt = 1, extent = [*xlims, *ylims], alpha = 0.5, zorder = 2)

    bins = np.arange(*xlims, 0.2)
    P16, P50, P84 = fplt.average_line(np.log10(L['FUV']), beta,bins)

    ax.fill_between(bins,P84,P16,color='k', alpha=0.15)
    ax.plot(bins,P50,ls='-',c='k', alpha=1.0, lw=1)

    # add observations

    ax.scatter(Bouwens['log10LFUV'], Bouwens['beta'], zorder = 4)


    ax.text(0.8, 0.9, f'$z={snaps[snap]}$', transform=ax.transAxes)




    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.set_ylabel(r'$\beta$')
    ax.set_xlabel(r'$\rm\log_{10}(L_{FUV}/erg\ s^{-1}\ Hz^{-1})$')
    ax.grid(True)
    fig.savefig('figs/beta_luminosity.pdf')






# --- colour-coded by attenuation

beta_luminosity_attenuation = False
if beta_luminosity_attenuation:

    xlims = [28, 30.5]
    ylims = [-2.75,-1.0]
    fig, ax, cax = fplt.single_wcbar()

    L = {}
    for f in ['FUV','NUV']:
        L[f] = np.array([])
        L[f'{f}_intrinsic'] = np.array([])
        for sim in data.keys():
            L[f] = np.hstack((L[f], data[f'{sim}/{snap}/Galaxy/BPASS/ChabrierIMF/Luminosity/DustModelI/{f}'][:]))
            L[f'{f}_intrinsic'] = np.hstack((L[f'{f}_intrinsic'], data[f'{sim}/{snap}/Galaxy/BPASS/ChabrierIMF/Luminosity/Intrinsic/{f}'][:]))


    beta = np.log10(L['FUV']/L['NUV'])/np.log10(1500./2500.) - 2.0

    AFUV = 2.5*np.log10(L['FUV_intrinsic']/L['FUV'])

    print(np.max(AFUV))

    norm = mpl.colors.Normalize(vmin=-3.0, vmax=3.0)

    ax.scatter(np.log10(L['FUV']), beta, c=cm.coolwarm(norm(AFUV)), s=1, alpha = 1.0, zorder = 1)

    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cm.coolwarm), cax=cax)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'$\rm A_{FUV}$', fontsize=8)

    ax.text(0.8, 0.9, f'$z={snaps[snap]}$', transform=ax.transAxes)

    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.set_ylabel(r'$\beta$')
    ax.set_xlabel(r'$\rm\log_{10}(L_{FUV}/erg\ s^{-1}\ Hz^{-1})$')
    ax.grid(True)
    fig.savefig('figs/beta_luminosity_attenuation.pdf')




# --- colour-coded by ssfr

beta_attenuation_ssfr = False
if beta_attenuation_ssfr:
    xlims = [-2.75, -1.0]
    ylims = [-0.5, 6.0]
    fig, ax, cax = fplt.single_wcbar()

    sfr = np.array([])
    mstar = np.array([])

    for sim in data.keys():
        sfr = np.hstack((sfr, data[f'{sim}/{snap}/Galaxy']['SFR']['SFR_10']))
        mstar = np.hstack((mstar, data[f'{sim}/{snap}/Galaxy']['Mstar_30']))


    s = np.log10(mstar)>8.5

    ssfr = np.log10(sfr)-np.log10(mstar)-1+9

    L = {}
    for f in ['FUV','NUV']:
        L[f] = np.array([])
        L[f'{f}_intrinsic'] = np.array([])
        for sim in data.keys():
            L[f] = np.hstack((L[f], data[f'{sim}/{snap}/Galaxy/BPASS/ChabrierIMF/Luminosity/DustModelI/{f}'][:]))
            L[f'{f}_intrinsic'] = np.hstack((L[f'{f}_intrinsic'], data[f'{sim}/{snap}/Galaxy/BPASS/ChabrierIMF/Luminosity/Intrinsic/{f}'][:]))


    beta = np.log10(L['FUV']/L['NUV'])/np.log10(1500./2500.) - 2.0

    AFUV = 2.5*np.log10(L['FUV_intrinsic']/L['FUV'])

    norm = mpl.colors.Normalize(vmin=-0.5, vmax=1.5)

    ax.scatter(beta[s], AFUV[s], c=cm.coolwarm(norm(ssfr[s])), s=1, alpha = 1.0, zorder = 1)


    # --- basic theory

    A1500 = np.arange(0,6,0.1)

    beta_int = -2.4

    beta = A1500*(1-(2500/1500)**-1) + beta_int

    print(beta)

    ax.plot(beta, A1500)


    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cm.coolwarm), cax=cax)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'$\log_{10}({\rm sSFR}/{\rm Gyr^{-1})}$', fontsize=8)

    ax.text(0.8, 0.9, f'$z={snaps[snap]}$', transform=ax.transAxes)

    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.set_xlabel(r'$\beta$')
    ax.set_ylabel(r'$\rm A_{FUV}$')
    ax.grid(True)
    fig.savefig('figs/beta_attenuation_ssfr.pdf')




# --- not color-coded

beta_attenuation = True
if beta_attenuation:
    xlims = [-2.75, -1.0]
    ylims = [-0.5, 6.0]
    fig, ax = fplt.single()


    sfr = np.array([])
    mstar = np.array([])

    for sim in data.keys():
        sfr = np.hstack((sfr, data[f'{sim}/{snap}/Galaxy']['SFR']['SFR_10']))
        mstar = np.hstack((mstar, data[f'{sim}/{snap}/Galaxy']['Mstar_30']))


    s = np.log10(mstar)>8.5



    L = {}
    for f in ['FUV','NUV']:
        L[f] = np.array([])
        L[f'{f}_intrinsic'] = np.array([])
        for sim in data.keys():
            L[f] = np.hstack((L[f], data[f'{sim}/{snap}/Galaxy/BPASS/ChabrierIMF/Luminosity/DustModelI/{f}'][:]))
            L[f'{f}_intrinsic'] = np.hstack((L[f'{f}_intrinsic'], data[f'{sim}/{snap}/Galaxy/BPASS/ChabrierIMF/Luminosity/Intrinsic/{f}'][:]))

    beta = np.log10(L['FUV']/L['NUV'])/np.log10(1500./2500.) - 2.0

    AFUV = 2.5*np.log10(L['FUV_intrinsic']/L['FUV'])

    ax.hexbin(beta[s], AFUV[s], gridsize = (25,25), bins = 'log', cmap='Greys', linewidths=0., mincnt = 1, extent = [*xlims, *ylims], alpha = 0.5, zorder = 2)


    # --- basic theory

    A1500 = np.arange(0,6,0.1)
    beta_int = -2.4
    beta = A1500*(1-(2500/1500)**-1) + beta_int
    ax.plot(beta, A1500)

    A1500 = np.arange(0,6,0.1)
    beta_int = -2.4
    beta = A1500*(1-(2500/1500)**-0.7) + beta_int
    ax.plot(beta, A1500)

    A1500 = np.arange(0,6,0.1)
    beta_int = -2.4
    beta = A1500*(1-(2500/1500)**-1.3) + beta_int
    ax.plot(beta, A1500)

    ax.text(0.8, 0.9, f'$z={snaps[snap]}$', transform=ax.transAxes)
    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.set_xlabel(r'$\beta$')
    ax.set_ylabel(r'$\rm A_{FUV}$')
    ax.grid(True)
    fig.savefig('figs/beta_attenuation.pdf')




# --- intrinsic vs. observed

comparison = False
if comparison:
    xlims = [28, 30.5]
    ylims = [-2.75,-1.0]
    fig, ax = fplt.single()

    StellarMass = np.array([])
    for sim in data.keys():
        StellarMass = np.hstack((StellarMass, data[f'{sim}/{snap}/Galaxy']['Mstar_30']))

    print(StellarMass.shape)

    s = np.log10(StellarMass) > 9.0

    L = {}
    for f in ['FUV','NUV']:
        L[f] = np.array([])
        L[f'{f}_intrinsic'] = np.array([])
        for sim in data.keys():
            L[f] = np.hstack((L[f], data[f'{sim}/{snap}/Galaxy/BPASS/ChabrierIMF/Luminosity/DustModelI/{f}'][:]))
            L[f'{f}_intrinsic'] = np.hstack((L[f'{f}_intrinsic'], data[f'{sim}/{snap}/Galaxy/BPASS/ChabrierIMF/Luminosity/Intrinsic/{f}'][:]))

    beta = np.log10(L['FUV']/L['NUV'])/np.log10(1500./2500.) - 2.0
    beta_intrinsic = np.log10(L['FUV_intrinsic']/L['NUV_intrinsic'])/np.log10(1500./2500.) - 2.0


    print(L['FUV'][s].shape)

    for l,l_int,b,b_int in zip(L['FUV'][s], L['FUV_intrinsic'][s], beta[s], beta_intrinsic[s]):

        ax.plot([np.log10(l_int), np.log10(l)],[b_int, b], c='k', alpha = 0.1)

    ax.text(0.8, 0.9, f'$z={snaps[snap]}$', transform=ax.transAxes)

    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.set_ylabel(r'$\beta$')
    ax.set_xlabel(r'$\rm\log_{10}(L_{FUV}/erg\ s^{-1}\ Hz^{-1})$')
    ax.grid(True)
    fig.savefig('figs/beta-comparison.pdf')
