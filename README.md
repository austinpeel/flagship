# flagship
Python codes to process and experiment with the `Euclid Flagship` mock using weak-lensing data.

Galaxy catalogs covering an octant of simulated sky can be downloaded from [CosmoHUB](https://cosmohub.pic.es/catalogs).  
Previous version: 1.3.3 (galaxies/lensing mismatch error)  
Current version:  1.5.2 (error fixed)

Although the official Flagship release is not available to non-Euclid members, some versions of the previous MICE catalog are (also found on CosmoHUB).

## Introduction
This package contains Python 2.7 modules I have written to access, process, and experiment with the weak-lensing mock data of the latest Euclid Flagship release. It was primarily created for personal use to help me test ideas related to estimating the masses of galaxy clusters using weak lensing. It is still a work in progress. As such, it is not (yet) documented as well as it could be, and I do not make any guarantees regarding its accuracy. However, if you happen to find any part of it useful in your own work, please feel free. I will also be happy to answer questions about how it works.

## Notes
The columns of the Flagship catalog used in this study are
1. `halo_id`
2. `galaxy_id`
3. `ra_gal`
4. `dec_gal`
5. `ra_gal_mag`
6. `dec_gal_mag`
7. `kappa`
8. `gamma1`
9. `gamma2`
10. `true_redshift_gal`
11. `observed_redshift_gal`
12. `kind`
13. `disk_angle`
14. `disk_axis_ratio`
15. `num_p`
16. `halo_lm`
17. `true_redshift_halo`

I have checked that indeed `disk_angle` = `bulge_angle` and
`disk_axis_ratio` = `bulge_axis_ratio`, so downloading only one of each
column is necessary for generating intrinsic galaxy ellipticities.

### Examples
TODO
```python
from wltools.mapping import bin2d
galmap = bin2d(cat.x_gal, cat.y_gal, npix=128)
```

### Questions
Why are the vast majority of `ra_gal_mag` and `dec_gal_mag` values NAN ?  
[Update] This has been fixed in the second release.

## Plans
The general plan is to test the weak-lensing masses of clusters in the simulation. This will likely be done by fitting NFW profiles to individual clusters, at least where they are massive enough (greater than a few x10^14 solar masses). Another possibility is to simply calculate aperture masses within a given radius.

### Step 1
Downloaded the mock catalog. To start with, I have selected a 10 x 10 deg^2
patch spanning RA, Dec in [0, 10) deg.

### Step 2
Verify that the DM halos trace the convergence field.
Based on tests so far (as of 4 Oct 2017), there seems to be an issue with this.
The convergence field does not track the galaxy distribution as it should.
With Arnau Pujol, we have calculated $\kappa_g$ from the weighted projection
of $\delta_g$ as in Eq. (16) of Pujol et al. (2016) [arXiv:1601.00160] using
all galaxies with z < 1.1. Comparing with the kappa map from galaxies in
1.0 < z < 1.2, there is essentially no correlation.

Pablo Fosalba, Jorge Carretero, and Romain Teyssier have been notified, and
we are all working to figure out what is wrong.

[Update] Version 1.5.2 released on Dec. 22, 2017 appears to have fixed this.
