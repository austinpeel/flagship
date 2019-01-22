import time
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from wltools.plotting import plot2d, plot_sticks
from wltools.mapping import bin2d, ks93
from wltools.projections import gnom
from wltools.utils import print_time
from wltools.conversions import to_deg
from .io import fetch_cat
from .octant import which_patch, which_patches, patch_extent


class flagship_field(object):
    def __init__(self, patch_id=None, extent=None, version='1.5.2', timed=True):
        start_time = time.time()
        if (patch_id is None) and (extent is None):
            print("Error: must provide either `patch` or `extent`.")
            return

        if patch_id is not None:
            extent = patch_extent(patch_id)
            ra0 = (extent[0] + extent[1]) / 2.
            dec0 = (extent[2] + extent[3]) / 2.
            patch_ids = [patch_id]
        elif extent is not None:
            extent = to_deg(extent)
            ra0 = (extent[0] + extent[1]) / 2.
            dec0 = (extent[2] + extent[3]) / 2.
            patch_ids = which_patches(extent)

        self.ra0 = ra0
        self.dec0 = dec0
        self.extent = extent
        self.aspect = ((extent[1] - extent[0]) / (extent[3] - extent[2])).value
        self.patch_ids = patch_ids

        # Load galaxies
        ra_gal = []
        dec_gal = []
        z_gal = []
        gamma1 = []
        gamma2 = []
        kappa = []
        e1_gal = []
        e2_gal = []
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

        # Flatten and store lists
        self.ra_gal = np.array([v for sublist in ra_gal for v in sublist])
        self.dec_gal = np.array([v for sublist in dec_gal for v in sublist])
        self.z_gal = np.array([v for sublist in z_gal for v in sublist])
        self.gamma1 = -np.array([v for sublist in gamma1 for v in sublist])
        self.gamma2 = -np.array([v for sublist in gamma2 for v in sublist])
        self.kappa = np.array([v for sublist in kappa for v in sublist])
        self.ngal = len(self.ra_gal)

        # Tangent plane projection
        self.x_gal, self.y_gal = gnom.radec2xy(self.ra0, self.dec0,
                                               self.ra_gal, self.dec_gal)
        # Geometry
        ll = gnom.radec2xy(self.ra0, self.dec0, self.extent[0], self.extent[2]) # lower left
        lc = gnom.radec2xy(self.ra0, self.dec0, self.ra0, self.extent[2])       # lower center
        lr = gnom.radec2xy(self.ra0, self.dec0, self.extent[1], self.extent[2]) # lower right
        ul = gnom.radec2xy(self.ra0, self.dec0, self.extent[0], self.extent[3]) # upper left
        ur = gnom.radec2xy(self.ra0, self.dec0, self.extent[1], self.extent[3]) # upper right
        delta_x = lr[0][0] - ll[0][0]
        delta_y = ul[1][0] - lc[1][0]
        self.aspect_xy = delta_x / delta_y
        self.extent_xy = [ll[0][0], lr[0][0], lc[1][0], ul[1][0]] * u.rad

        self.version = version

        if timed:
            printtime(time.time() - start_time)


    def plot_kappa(self, npix=128, ks=False, proj=False, lm_min=None,
                   zmax=None, **kwargs):
        if proj:
            x_gal = self.x_gal
            y_gal = self.y_gal
            extent = self.extent_xy.value
            unit = self.extent_xy.unit
        else:
            x_gal = self.ra_gal
            y_gal = self.dec_gal
            extent = self.extent.value
            unit = self.extent.unit

        gamma1 = self.gamma1
        gamma2 = self.gamma2
        kappa = self.kappa

        if zmax is not None:
            sel = self.z_gal < zmax
            x_gal = x_gal[sel]
            y_gal = y_gal[sel]
            gamma1 = gamma1[sel]
            gamma2 = gamma2[sel]
            kappa = kappa[sel]

        if ks:
            g1map, g2map = bin2d(x_gal, y_gal, v=[gamma1, gamma2],
                                 npix=npix, extent=extent)
            kappaE, kappaB = ks93(g1map, g2map)
        else:
            kappaE = bin2d(x_gal, y_gal, v=kappa, npix=npix,
                           extent=extent)
            kappaB = np.zeros_like(kappaE)

        print("< kappa > = {:.4e}".format(kappaE.mean()))

        fig, ax = plt.subplots(1, 1, figsize=(6.5, 6))
        plot2d(kappaE, fig=fig, ax=ax, extent=extent, **kwargs)
        ax.set_xlabel("RA [{}]".format(unit))
        ax.set_ylabel("Dec [{}]".format(unit))
        ax.set_title("{} $\kappa$".format(["true", "KS93"][ks]))

        if lm_min is not None:
            for pid in self.patch_ids:
                halos = fetch_cat(pid, halos=True, version=self.version)
                if halos is None:
                    continue
                sel = ((halos['lmfof'] > lm_min) &
                       (halos['ra'] > extent[0]) &
                       (halos['ra'] < extent[1]) &
                       (halos['dec'] > extent[2]) &
                       (halos['dec'] < extent[3]))
                for halo in halos[sel]:
                    x = halo['ra']
                    y = halo['dec']
                    z = halo['z']
                    ax.scatter(x, y, s=20, c='w', alpha=0.6)
                    dy = 0.02 * (extent[3] - extent[2])
                    ax.text(x, y + dy, "{:.2f}".format(z),
                            color='w', fontsize=6, ha='center')

        if ks:
            fig, ax = plt.subplots(1, 1, figsize=(6.5, 6))
            plot2d(kappaB, fig=fig, ax=ax, extent=extent, **kwargs)
            ax.set_xlabel("RA [{}]".format(unit))
            ax.set_ylabel("Dec [{}]".format(unit))
            ax.set_title("KS93 $\kappa_\mathrm{B}$")


    def plot_halos(self, lm_min=14.0):
        fig, ax = plt.subplots(1, 1, figsize=(6.5, 6))
        for pid in self.patch_ids:
            halos = fetch_cat(patch_id=pid, halos=True, version=self.version)
            if halos is None:
                continue
            extent = self.extent.value
            sel = ((halos['lmfof'] > lm_min) &
                   (halos['ra'] > extent[0]) &
                   (halos['ra'] < extent[1]) &
                   (halos['dec'] > extent[2]) &
                   (halos['dec'] < extent[3]))
            for halo in halos[sel]:
                ax.scatter(halo['ra'], halo['dec'], s=5, c='k', alpha=0.6)

        ax.set_xlabel('RA')
        ax.set_ylabel("Dec")
        ax.set_xlim(extent[:2])
        ax.set_ylim(extent[2:])
