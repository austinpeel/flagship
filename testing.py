import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
from scipy.signal import correlate
from wltools.plotting import plot2d
from wltools.mapping import bin2d
from wltools.projections import gnom
from wltools.stats import powspec
from .io import *
from .utils import split_extent


def kappa_x_halos(split_id=0, npix=512, save=False):
    """
    Check that the halo distribution traces (at least roughly) the
    convergence map.
    """
    halofile = 'splits/{}/halos/halocat_13.5.fits'.format(split_id)
    halos = fetch_cat(halofile)

    # Restrict halos
    halo_zmin = 0
    halo_zmax = 1
    halo_lmmin = 14
    sel = ((halos['z'] > halo_zmin) & (halos['z'] < halo_zmax) &
           (halos['lmvir'] > halo_lmmin))
    halos = halos[sel]
    print("Halos")
    print("-----")
    print("{} with {:.2f} < z < {:.2f}".format(sel.sum(), halo_zmin, halo_zmax))
    print("    log10(m_vir) > {:.1f}".format(halo_lmmin))
    ra_halo = halos['ra'] * u.deg
    dec_halo = halos['dec'] * u.deg

    # Galaxies
    cat = fetch_cat('splits/{}/galcat_full.fits'.format(split_id))
    # cat = cat[cat['true_redshift_gal'] < 1.5]
    ra_gal = cat['ra_gal_mag'] * u.deg
    dec_gal = cat['dec_gal_mag'] * u.deg

    # Projections
    ramin, ramax, decmin, decmax = split_extent(split_id)
    ra0 = 0.5 * (ramin + ramax) # [deg]
    dec0 = 0.5 * (decmin + decmax) # [deg]
    x_halo, y_halo = gnom.radec2xy(ra0, dec0, ra_halo, dec_halo) * u.rad
    x_gal, y_gal = gnom.radec2xy(ra0, dec0, ra_gal, dec_gal) * u.rad

    x_halo = ra0 + x_halo.to(u.deg)
    y_halo = dec0 + y_halo.to(u.deg)
    x_gal = ra0 + x_gal.to(u.deg)
    y_gal = dec0 + y_gal.to(u.deg)

    xmin, y = gnom.radec2xy(ra0, dec0, ramin, decmin) * u.rad
    xmax, y = gnom.radec2xy(ra0, dec0, ramax, decmin) * u.rad
    x, ymin = gnom.radec2xy(ra0, dec0,  ra0, decmin) * u.rad
    x, ymax = gnom.radec2xy(ra0, dec0, ramin, decmax) * u.rad

    xmin = ra0 + xmin.to(u.deg)
    xmax = ra0 + xmax.to(u.deg)
    ymin = dec0 + ymin.to(u.deg)
    ymax = dec0 + ymax.to(u.deg)

    fig, ax = plt.subplots(1, 1, facecolor='w')
    size = 60 * (halos['lmvir'] - halo_lmmin) + 1
    ax.scatter(x_halo.value, y_halo.value, facecolors='none', edgecolors='r',
               s=size, alpha=0.6)
    ax.set_aspect('equal')
    ax.set_xlim(ramin.value, ramax.value)
    ax.set_ylim(decmin.value, decmax.value)
    ax.set_xlabel('RA [deg]')
    ax.set_ylabel('Dec [deg]')

    ext = ax.get_xlim() + ax.get_ylim()
    kappa = bin2d(x_gal.value, y_gal.value, v=cat['kappa'], npix=npix)
    plot2d(gaussian_filter(kappa, 0), fig=fig, ax=ax, cmap='bone',
           extent=ext)

    # ax.set_xlim(0, xmax.value)
    # ax.set_ylim(ymin.value, 0)
    ax.set_xlim(ra0.value, ramax.value)
    ax.set_ylim(decmin.value, dec0.value)
    # ax.set_xlim(ramin.value, ra0.value)
    # ax.set_ylim(decmin.value, dec0.value)

    plt.show()

    if save:
        outfile = datapath('splits/{}/halos_x_kappa.png'.format(split_id))
        plt.savefig(outfile, format='png')


def xcorr(split_id=0, npix=512):
    cat = fetch_cat('splits/{}/galcat_full.fits'.format(split_id))
    # Source plane
    sel_src = ((cat['true_redshift_gal'] > 0.95) &
               (cat['true_redshift_gal'] < 1.05))
    ra_src = cat['ra_gal_mag'][sel_src]
    dec_src = cat['dec_gal_mag'][sel_src]
    kappa_src = cat['kappa'][sel_src]

    kappa = bin2d(ra_src, dec_src, v=kappa_src, npix=npix)
    # plot2d(kappa, cmap='magma')
    print("<kappa> = {:.3f}".format(kappa.mean()))

    # Lens plane
    halofile = 'splits/{}/halos/halocat_13.5.fits'.format(split_id)
    halos = fetch_cat(halofile)
    lm_min = 13.8
    sel_halos = ((halos['z'] > 0.45) & (halos['z'] < 0.55) &
                 (halos['lmvir'] > lm_min))
    halos = halos[sel_halos]
    print("{} with {:.2f} < z < {:.2f}".format(sel_halos.sum(), 0.45, 0.55))
    print("    log10(m_vir) > {:.1f}".format(lm_min))
    ra_halo = halos['ra']
    dec_halo = halos['dec']
    fig, ax = plt.subplots(1, 1, facecolor='w')
    size = 60 * (halos['lmvir'] - lm_min) + 1
    ax.scatter(ra_halo, dec_halo, facecolors='none', edgecolors='r',
               s=size, alpha=0.6)
    ax.set_aspect('equal')
    ramin, ramax, decmin, decmax = split_extent(split_id).value
    ax.set_xlim(ramin, ramax)
    ax.set_ylim(decmin, decmax)

    sel_lens = ((cat['true_redshift_gal'] > 0.45) &
                (cat['true_redshift_gal'] < 0.55))
    ra_lens = cat['ra_gal_mag'][sel_lens]
    dec_lens = cat['dec_gal_mag'][sel_lens]
    delta_lens = bin2d(ra_lens, dec_lens, npix=npix)
    delta = delta_lens - delta_lens.mean()
    print("<delta> = {:.3f}".format(delta.mean()))
    plot2d(delta, fig=fig, ax=ax, cmap='bone', extent=(ramin, ramax, decmin, decmax))

    # corr = correlate(kappa, kappa, method='fft', mode='same')
    l, cl = powspec(kappa * delta, 10)
    # plot2d(corr)
    fig2, ax2 = plt.subplots(1, 1)
    ax2.loglog(l, l * (l + 1) * cl / (2 * np.pi))


def test_mice(save=False):
    from astropy.io import fits

    ramin = 0
    ramax = 10
    decmin = 60
    decmax = 65
    # --------
    #  Halos
    # --------
    halos = fits.getdata('/Users/apeel/Data/CFC4/M2_mice2_500deg2mice2_500_no_masked_halo.fits')
    lm_min = 14
    sel = ((halos['ra'] > ramin) & (halos['ra'] < ramax) &
           (halos['dec'] > decmin) & (halos['dec'] < decmax) &
        #    (halos['z_cos'] > 0.35) & (halos['z_cos'] < 0.45) &
           (halos['mhhalo'] > lm_min))
    halos = halos[sel]
    ra_halo = halos['ra'] * u.deg
    dec_halo = halos['dec'] * u.deg

    # -----------
    #  Galaxies
    # -----------
    cat = fits.getdata('/Users/apeel/Data/CFC4/splits/join_0_1_13_14.fits')
    sel = (cat['dec_gal'] < decmax)
    # sel = sel & (cat['z_gal'] > 0.45) & (cat['z_gal'] < 1.05)
    cat = cat[sel]
    ra_gal = cat['ra_gal'] * u.deg
    dec_gal = cat['dec_gal'] * u.deg

    # ------------
    #  Projection
    # ------------
    ra0 = np.mean([ramin, ramax]) * u.deg
    dec0 = np.mean([decmin, decmax]) * u.deg
    x_halo, y_halo = gnom.radec2xy(ra0, dec0, ra_halo, dec_halo) * u.rad
    x_gal, y_gal = gnom.radec2xy(ra0, dec0, ra_gal, dec_gal) * u.rad

    # x_halo = x_halo.to(u.arcmin).value
    # y_halo = y_halo.to(u.arcmin).value
    # x_gal = x_gal.to(u.arcmin).value
    # y_gal = y_gal.to(u.arcmin).value

    xmin, y = gnom.radec2xy(ra0, dec0, ramin * u.deg, decmin * u.deg)
    xmax, y = gnom.radec2xy(ra0, dec0, ramax * u.deg, decmin * u.deg)
    x, ymin = gnom.radec2xy(ra0, dec0,  ra0, decmin * u.deg)
    x, ymax = gnom.radec2xy(ra0, dec0, ramin, decmax * u.deg)
    # ext = [xmin[0], xmax[0], ymin[0], ymax[0]]

    # ------
    #  Plot
    # ------
    fig, ax = plt.subplots(1, 1, facecolor='w')
    size = 60 * (halos['mhhalo'] - lm_min) + 1
    ax.scatter(x_halo, y_halo, facecolors='none', edgecolors='r',
               s=size, alpha=0.6)
    ax.scatter(0, 0, marker='x', c='w', alpha=0.6)
    ax.set_aspect('equal')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel('$\Delta$RA [rad]')
    ax.set_ylabel('$\Delta$Dec [rad]')
    ax.set_title("MICE CFC4 $\kappa$ patch")

    ext = ax.get_xlim() + ax.get_ylim()
    kappa = bin2d(x_gal.value, y_gal.value, v=cat['kappa'], npix=256)
    plot2d(gaussian_filter(kappa, 0), fig=fig, ax=ax, cmap='bone',
           extent=ext)

    if save:
        outfile = '/Users/apeel/Desktop/cfc4_halos_x_kappa.png'
        plt.savefig(outfile, format='png')
