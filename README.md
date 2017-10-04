# flagship
Python codes to process the `Euclid Flagship` mock and to do weak-lensing analyses.

The catalog can be downloaded from [CosmoHUB](https://cosmohub.pic.es/#/catalogs/53)
---

## Introduction
Testing...

```python
from wltools.mapping import bin2d
galmap = bin2d(cat.x_gal, cat.y_gal, npix=128)
```

---

## Notes
Columns used in this study are
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
12. `disk_angle`
13. `disk_axis_ratio`
14. `num_p`
15. `halo_lm`
16. `true_redshift_halo`

I have checked that indeed `disk_angle` = `bulge_angle` and
`disk_axis_ratio` = `bulge_axis_ratio`, so downloading only one of each
column is necessary for generating intrinsic galaxy ellipticities.

### Questions
Why are the vast majority of `ra_gal_mag` and `dec_gal_mag` values NAN ?
[Update] This has been fixed in the second release.

---

## Plans
The plan is to test the `weak-lensing masses` of clusters in the simulation.

### Step 1
Download the mock catalog. To start with, I have selected a 10 x 10 deg^2
patch spanning RA, Dec in [0, 10) deg.

### Step 2
Verify that the DM halos trace the convergence field.
