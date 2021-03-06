import os
import time
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from scipy.ndimage.filters import gaussian_filter
from shapely.geometry import Point
from wltools.conversions import to_arcmin
from wltools.projections import gnom
from wltools.mapping import bin2d
from wltools.plotting import plot2d, plot_sticks
from wltools.halos.fitting import lsqfit
from wltools.utils import printtime
from .io import fetch_cat
from .field import flagship_field, which_patch, which_patches
from .cosmo import flagship_cosmology


def rank2id(rank, patch_id, version='1.5.2'):
    """
    Return the unique halo_id of the `n`th most massive halo by fof mass,
    starting with 0. Only halos with log(m_fof) > 13.5 are considered.
    """
    halos = fetch_cat(patch_id, halos=True, version=version)
    nhalos = len(halos)
    if rank > nhalos:
        print("There are only {} halos with log(m_fof) > 13.5.".format(nhalos))
        halo_id = -1
    else:
        halo_id = halos[::-1]['halo_id'][rank]
    return halo_id


class flagship_halo(object):
    """Euclid Flagship dark matter halo class."""
    def __init__(self, patch_id, halo_id=None, rank=0, version='1.5.2'):
        """
        Parameters
        ----------
        patch_id : int
            ID of the patch (0-80) the halo lies in. This is necessary, since
            we generate halo catalogues from each patch of the galaxy
            catalogue; i.e. we need to know where the halo lives.
        halo_id : int, optional
            Unique halo identifier given in Flagship catalogue.
            Default is None.
        rank : int, optional
            Ranking by log(m_fof) of halos in area corresponding to patch_id.
            Default is 0, i.e. the most massive halo.

        Note
        ----
        One of either halo_id or rank must be provided, with halo_id
        taking precedence.
        """
        if halo_id is None:
            halo_id = rank2id(rank=rank, patch_id=patch_id, version=version)
            if halo_id < 0:
                print("Error: invalid rank.")
                return

        halos = fetch_cat(patch_id=patch_id, halos=True, version=version)
        if halo_id not in halos['halo_id']:
            print("Invalid halo id.")
            return

        # Extract the correct halo
        halo = halos[halos['halo_id'] == halo_id]

        # Cosmology
        self.cosmo = flagship_cosmology()

        # Halo properties (with units where appropriate)
        self.patch_id = patch_id
        self.id = halo['halo_id'][0]
        self.ra = halo['ra'][0] * u.degree
        self.dec = halo['dec'][0] * u.degree
        self.z = halo['z'][0]
        self.lmfof = halo['lmfof'][0] # [log10(M_sun / h)]
        self.mfof = 10**self.lmfof / self.cosmo.h * u.solMass
        if 'lmvir' in halos.names:
            self.lmvir = halo['lmvir'][0] # [log10(M_sun / h)]
            self.mvir = 10**self.lmvir / self.cosmo.h * u.solMass
        if 'lm200c' in halos.names:
            self.lm200c = halo['lm200c'][0] # [log10(M_sun / h)]
            self.m200c = 10**self.lm200c / self.cosmo.h * u.solMass
        self.ngal = halo['ngal'][0]
        self.npart = halos['npart'][0]
        self.version = version

    def load_members(self):
        # Load full galaxy catalog
        cat = fetch_cat(self.patch_id, version=self.version)
        sel = cat['halo_id'] == self.id
        print(sel.sum(), self.ngal)
        # if not sel.sum() == self.ngal:
            # print("Not enough. Need to load nearby patches.")
        print("THIS CODE IS UNFINISHED")
        return
        # Load galaxies
        ra_gal = []
        dec_gal = []
        z_gal = []
        for pid in patch_ids:
            print("Fetching patch {}".format(pid))
            cat = fetch_cat(pid, version=version)
            if cat is None:
                print("!! Skipped patch {} !! Need to download.".format(pid))
                continue
            print("Selecting galaxies.")
            # selection_time = time.time()
            sel = ((cat['ra_gal_mag'] > self.extent[0].value) &
                   (cat['ra_gal_mag'] < self.extent[1].value) &
                   (cat['dec_gal_mag'] > self.extent[2].value) &
                   (cat['dec_gal_mag'] < self.extent[3].value))
            print("Loaded {} galaxies within extent.".format(sel.sum()))
            # printtime(time.time() - selection_time)
            cat = cat[sel]
            ra_gal.append(cat['ra_gal_mag']) # a bit faster than extend()
            dec_gal.append(cat['dec_gal_mag'])
            z_gal.append(cat['true_redshift_gal'])
            gamma1.append(cat['gamma1'])
            gamma2.append(cat['gamma2'])
            kappa.append(cat['kappa'])
        return
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
    """A lensing postage stamp centered on a Flagship halo."""
    def __init__(self, patch_id, halo_id=None, rank=0, size=None,
                 version='1.5.2', timed=True):
        """
        Parameters
        ----------
        patch_id : int
            ID of the patch (0-80) the halo lives in.
        halo_id : int, optional
            Unique halo id.
        rank : int, optional
            Ranking by log(m_fof) of halos in area corresponding to patch_id.
            Default is 0, i.e. the most massive halo.
        size : float, optional (as `astropy.units.quantity.Quantity`)
            Angular size of square postage stamp in arcmin, centered on the
            halo. If no compatible units are given, arcmin are assumed.
            Default is None, in which case the size is computed as the
            angular extent corresponding to 8 Mpc at the halo redshift.
        """
        starttime = time.time()

        # Load halo properties as a flagship_halo object
        self.halo = flagship_halo(patch_id=patch_id, halo_id=halo_id,
                                  rank=rank, version=version)
        self.patch_id = patch_id
        self.version = version

        # Cosmology
        self.cosmo = self.halo.cosmo
        angular_scale = self.cosmo.arcsec_per_kpc_proper(self.halo.z)
        self.angular_scale = angular_scale.to(u.arcmin / u.Mpc)

        # Set stamp size (in radians)
        if size is not None:
            try:
                size = to_arcmin(size).to('rad')
            except u.UnitConversionError:
                print("Error: Invalid unit for size. Resorting to Default.")
                size = (8 * u.Mpc) * self.angular_scale
                size = size.to(u.rad)

        else:
            size = (8 * u.Mpc) * self.angular_scale
            size = size.to(u.rad)

        # Determine minimum rectangle in RA/Dec to contain the projected patch
        radec_extent = gnom.tangent_square(self.halo.ra, self.halo.dec, size)

        #
        ff = flagship_field(extent=radec_extent)

        # Pull all galaxies within the postage stamp area
        cat = fetch_cat(self.halo.patch_id)
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
        self._load_nearby_halos()

        if timed:
            printtime(time.time() - starttime)


    def _load_nearby_halos(self):
        halos = fetch_cat(self.patch_id, halos=True, version=self.version)
        # Project halos and select the ones within this halo's extent
        x_halo, y_halo = gnom.radec2xy(self.halo.ra, self.halo.dec,
                                       halos['ra'] * u.deg,
                                       halos['dec'] * u.deg)
        x_halo = (x_halo * u.rad).to(u.arcmin)
        y_halo = (y_halo * u.rad).to(u.arcmin)
        xmin, xmax, ymin, ymax = self.extent
        sel = ((x_halo > xmin) & (x_halo < xmax) &
               (y_halo > ymin) & (y_halo < ymax))
        halo_ids = halos['halo_id'][sel]
        self.nearby_halos = [flagship_halo(self.patch_id, hid) for hid
                             in halo_ids if hid != self.halo.id]


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


    def plot_nearby_halos(self, lm_min=None, ax=None, show_z=False):
        if ax is None:
            fig, ax = plt.subplots(1, 1, facecolor='w')
            ax.scatter(0, 0, marker='x')
        if lm_min is None:
            lm_min = 0.
        for halo in self.nearby_halos:
            if halo.lmfof > lm_min:
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
