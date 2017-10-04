import os
import time
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from scipy.ndimage.filters import gaussian_filter
from shapely.geometry import Point
from wltools.projections import gnom
from wltools.mapping import bin2d
from wltools.plotting import plot2d, plot_sticks
from wltools.halos.fitting import lsqfit
from .io import fetch_cat


class flagship_halo(object):
    def __init__(self, halo_id, split_id=0):
        # Cosmology
        Om0 = 0.319
        H0 = 67. # [km/s/Mpc]
        self.cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

        halofile = 'splits/{}/halos/halocat_13.5.fits'.format(split_id)
        halos = fetch_cat(halofile)
        # Halo properties (with units where appropriate)
        if halo_id not in halos['halo_id']:
            print("Invalid halo id.")
            return

        # Extract the correct halo
        halo = halos[halos['halo_id'] == halo_id]

        self.split_id = split_id
        self.id = halo['halo_id'][0]
        self.ra = halo['ra'][0] * u.degree
        self.dec = halo['dec'][0] * u.degree
        self.z = halo['z'][0]
        self.lmfof = halo['lmfof'][0] # [log10(M_sun / h)]
        self.lmvir = halo['lmvir'][0] # [log10(M_sun / h)]
        self.lm200c = halo['lm200c'][0] # [log10(M_sun / h)]
        self.mfof = 10**self.lmfof / self.cosmo.h  * u.solMass
        self.mvir = 10**self.lmvir / self.cosmo.h  * u.solMass
        self.m200c = 10**self.lm200c / self.cosmo.h * u.solMass
        self.ngal = halo['ngal'][0]
        self.npart = halos['npart'][0]

    def load_members(self):
        # Load full galaxy catalog
        cat = fetch_cat('splits/{}/galcat_full.fits'.format(self.split_id))
        # Locate entries matching unique halo index
        self.member_inds = (cat['halo_id'] == self.id)
        cat = cat[self.member_inds]
        self.id_gal = cat['galaxy_id']
        self.ra_gal = cat['ra_gal_mag'] * u.degree
        self.dec_gal = cat['dec_gal_mag'] * u.degree
        self.z_gal = cat['true_redshift_gal']
        # Project to tangent plane
        x_gal, y_gal = gnom.radec2xy(self.ra, self.dec,
                                     self.ra_gal, self.dec_gal)
        self.x_gal = (x_gal * u.rad).to(u.arcmin)
        self.y_gal = (y_gal * u.rad).to(u.arcmin)

    def projected_position(self, ra0, dec0):
        x, y = gnom.radec2xy(ra0, dec0, self.ra, self.dec)
        return ((x[0], y[0]) * u.rad).to(u.arcmin)


class postage_stamp(object):
    """Extract a lensing postage stamp centered on a Flagship halo."""
    def __init__(self, halo_id, size=None, timed=True):
        """
        Parameters
        ----------
        halo_id : int
            Unique halo id.
        size : float, optional (as `astropy.units.quantity.Quantity`)
            Angular size of square postage stamp in arcmin, centered on the
            halo. If no compatible units are given, arcmin are assumed.
            Default is None, in which case the size is computed as the
            angular extent corresponding to 8 Mpc at the halo redshift.
        """
        starttime = time.time()

        # Load halo properties as a flagship_halo object
        self.halo = flagship_halo(halo_id)

        # Cosmology
        self.cosmo = self.halo.cosmo
        angular_scale = self.cosmo.arcsec_per_kpc_proper(self.halo.z)
        self.angular_scale = angular_scale.to(u.arcmin / u.Mpc)

        # TODO This check can surely be improved...
        if size is not None:
            try:
                size = size.to(u.rad)
            except AttributeError:
                size = float(size) * u.arcmin
                size = size.to(u.rad)
            except u.UnitConversionError:
                print("Error: Invalid unit for size. Resorting to Default.")
                size = (8 * u.Mpc) * self.angular_scale
                size = size.to(u.rad)
        else:
            size = (8 * u.Mpc) * self.angular_scale
            size = size.to(u.rad)

        # Determine minimum rectangle in RA/Dec to contain the projected patch
        radec_extent = gnom.tangent_square(self.halo.ra, self.halo.dec, size)

        # Pull all galaxies within the postage stamp area
        cat = fetch_cat('splits/0/galcat_full.fits')
        halfsize = size.value / 2 # [radians]
        x, y = gnom.radec2xy(self.halo.ra, self.halo.dec,
                             cat['ra_gal_mag'], cat['dec_gal_mag'])
        sel = ((x >= -halfsize) & (x <= halfsize) &
               (y >= -halfsize) & (y <= halfsize))

        # Make catalog selection
        self.halo.load_members()
        self.member_mask = (sel & self.halo.member_inds)[sel]
        cat = cat[sel]
        self.id_gal = cat['galaxy_id']
        self.ra_gal = cat['ra_gal_mag'] * u.degree
        self.dec_gal = cat['dec_gal_mag'] * u.degree
        self.x_gal = (x[sel] * u.rad).to(u.arcmin)
        self.y_gal = (y[sel] * u.rad).to(u.arcmin)
        self.kappa = cat['kappa']
        self.gamma1 = -cat['gamma1']
        self.gamma2 = -cat['gamma2']
        self.z_gal = cat['true_redshift_gal']
        self.ngal = len(self.id_gal)
        self.area = (size.to(u.arcmin)**2)
        self.gal_density = self.ngal / self.area # arcmin^(-2)

        # Field properties
        radius = size.to(u.arcmin).value / 2
        self.extent = [-radius, radius, -radius, radius] * u.arcmin

        # Retrieve indices for background sources, i.e. galaxies that are
        # behind the lens and not members of the halo
        dz = 0
        self.src_mask = (self.z_gal > self.halo.z + dz) & (~self.member_mask)

        # Compute source density
        self.src_density = self.src_mask.sum() / self.area

        # Load nearby halos
        self.load_nearby_halos()

        if timed:
            print("Time : {0:.2f} sec".format(time.time() - starttime))


    def load_nearby_halos(self):
        halos = fetch_cat('splits/split0_halos_13.5.fits')
        # Project halos and select the ones within this halo's extent
        x_halo, y_halo = gnom.radec2xy(self.halo.ra, self.halo.dec,
                                       halos['halo_ra'] * u.deg,
                                       halos['halo_dec'] * u.deg)
        x_halo = (x_halo * u.rad).to(u.arcmin)
        y_halo = (y_halo * u.rad).to(u.arcmin)
        xmin, xmax, ymin, ymax = self.extent
        sel = ((x_halo > xmin) & (x_halo < xmax) &
               (y_halo > ymin) & (y_halo < ymax))
        halo_ids = halos['halo_id'][sel]
        self.nearby_halos = [flagship_halo(idx) for idx in halo_ids if
                             idx != self.halo.id]


    def plot_halo_members(self):
        print("{} member galaxies".format(self.halo.ngal))
        fig, ax = plt.subplots(1, 1, facecolor='w')
        sel = self.member_mask
        ax.scatter(self.halo.x_gal, self.halo.y_gal, s=2)
        ax.scatter(self.x_gal[sel], self.y_gal[sel], s=3)
        ax.scatter(0, 0, marker='x')
        ax.set_xlabel("$\Delta$RA [arcmin]")
        ax.set_ylabel("$\Delta$Dec [arcmin]")
        ax.set_xlim(self.extent[:2].value)
        ax.set_ylim(self.extent[-2:].value)
        ax.set_aspect('equal')
        plt.show()


    def plot_nearby_halos(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots(1, 1, facecolor='w')
            ax.scatter(0, 0, marker='x')
        for halo in self.nearby_halos:
            x, y = halo.projected_position(self.halo.ra, self.halo.dec)
            ax.scatter(x, y, c='r')
        ax.set_xlim(self.extent[:2].value)
        ax.set_ylim(self.extent[-2:].value)
        ax.set_aspect('equal')


    def _select_galaxies(self, foreground):
        if not foreground:
            sel = self.src_mask
        else:
            sel = np.ones(self.ngal) == np.ones(self.ngal)
        return sel, self.x_gal[sel].value, self.y_gal[sel].value


    def plot_density(self, npix=64, fg=False, save=False):
        sel, x_gal, y_gal = self._select_galaxies(foreground=fg)
        galmap = bin2d(x_gal, y_gal, npix=npix)
        print("< gal / pix > = {:.2f}".format(galmap.mean()))
        fig, ax = plt.subplots(1, 1, facecolor='w')
        plot2d(galmap, fig=fig, ax=ax, cmap='bone', extent=self.extent.value)
        ax.scatter(0, 0, marker='x', color='r', s=30, lw=1.5)
        ax.set_aspect('equal')
        ax.set_xlabel("$\Delta$RA [arcmin]")
        ax.set_ylabel("$\Delta$Dec [arcmin]")
        ax.set_title("{} galaxy density".format(['source', 'full'][fg]))
        plt.show()


    def plot_kappa(self, npix=64, fg=False, ker=0, save=False):
        sel, x_gal, y_gal = self._select_galaxies(foreground=fg)
        kappa = self.kappa[sel]
        galmap = bin2d(x_gal, y_gal, npix=npix)
        print("< gal / pix > = {:.2f}".format(galmap.mean()))
        kappamap = bin2d(x_gal, y_gal, v=kappa, npix=npix)
        print("< kappa > = {:.4f}".format(kappamap.mean()))
        fig, ax = plt.subplots(1, 1, facecolor='w')
        plot2d(gaussian_filter(kappamap, ker), fig=fig, ax=ax,
                 extent=self.extent.value, cmap='magma')
        ax.set_aspect('equal')
        ax.set_xlabel("$\Delta$RA [arcmin]")
        ax.set_ylabel("$\Delta$Dec [arcmin]")
        ax.set_title("{} convergence $\kappa$".format(['source', 'full'][fg]))
        plt.show()


    def plot_gamma(self, binned=True, npix=64, fg=False, scl=0.1, save=False):
        sel, x_gal, y_gal = self._select_galaxies(foreground=fg)
        masked = True
        if masked:
            nmax = int(30 * self.area.value)
            mask = np.random.choice(sel.sum(), nmax, replace=False)
            print("density = {:.2f}".format(nmax / self.area))
        else:
            mask = range(sel.sum())
        galmap = bin2d(x_gal[mask], y_gal[mask], npix=npix)
        print("< gal / pix > = {:.2f}".format(galmap.mean()))
        extent = self.extent.value
        if binned:
            gammas = [self.gamma1[sel][mask], self.gamma2[sel][mask]]
            g1map, g2map = bin2d(x_gal[mask], y_gal[mask], v=gammas, npix=npix)
            xx = np.linspace(extent[0], extent[1], npix)
            yy = np.linspace(extent[2], extent[3], npix)
            xv, yv = np.meshgrid(xx, yy, sparse=False)
            x = xv.flatten()
            y = yv.flatten()
            g1 = g1map.flatten()
            g2 = g2map.flatten()
        else:
            x = x_gal[mask]
            y = y_gal[mask]
            g1 = self.gamma1[sel][mask]
            g2 = self.gamma2[sel][mask]
        fig, ax = plt.subplots(1, 1, facecolor='w')
        plot_sticks(x, y, g1, g2, scl=scl, ax=ax, color='k', extent=extent)
        ax.scatter(0, 0, marker='x')
        rmin = ((750 * u.kpc) * self.angular_scale).to(u.arcmin).value
        rmax = ((3.0 * u.Mpc) * self.angular_scale).to(u.arcmin).value
        circle_min = Point(0, 0).buffer(rmin)
        circle_max = Point(0, 0).buffer(rmax)
        xc_min, yc_min = circle_min.exterior.xy
        xc_max, yc_max = circle_max.exterior.xy
        ax.plot(xc_min, yc_min, c='m', linewidth=1)
        ax.plot(xc_max, yc_max, c='m', linewidth=1)
        ax.set_aspect('equal')
        ax.set_xlabel("$\Delta$RA [arcmin]")
        ax.set_ylabel("$\Delta$Dec [arcmin]")
        ax.set_title("{} shear $\gamma$".format(['source', 'full'][fg]))

        self.plot_nearby_halos(ax=ax)
        plt.show()


    def fit(self, prof='nfw'):
        print("Performing least squares fit to {} profile.".format(prof))
        fit = lsqfit(self.halo.z)
        sel, x_gal, y_gal = self._select_galaxies(foreground=False)
        nmax = int(30 * self.area.value)
        mask = np.random.choice(sel.sum(), nmax, replace=False)
        x_gal = x_gal * u.arcmin
        y_gal = y_gal * u.arcmin
        e1_gal = self.gamma1[sel]
        e2_gal = self.gamma2[sel]
        z_gal = self.z_gal[sel]
        fit.load_galaxies(x_gal[mask], y_gal[mask], e1_gal[mask],
                          e2_gal[mask], z_gal[mask])
        return fit.do_fit()
