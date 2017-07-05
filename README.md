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
5. `kappa`
6. `gamma1`
7. `gamma2`
8. `true_redshift_gal`
9. `disk_angle`
10. `disk_axis_ratio`
11. `halo_lm`
12. `true_redshift_halo`

I have checked that indeed `disk_angle` = `bulge_angle` and
`disk_axis_ratio` = `bulge_axis_ratio`, so downloading only one of each
column is necessary for generating intrinsic galaxy ellipticities.

### Questions
Why are the vast majority of `ra_gal_mag` and `dec_gal_mag` values NAN ??

---

## Plans
The plan is to test the `weak-lensing masses` of clusters in the simulation.

### Step 1
Download the mock catalog.

### Step 2
Measure masses. Easy.
