import numpy as np
import astropy.units as u
from wltools.conversions import to_deg


def which_patch(ra, dec):
    """
    Determine which patch of the octant (numbered 0-80) the point (ra, dec)
    lies in. Patches are 100 deg^2 squares with boundaries defined as
            RA        Dec
    0  -- [0, 10)   [0, 10)
    1  -- [10, 20)  [0, 10)
           ...       ...
    8  -- [80, 90)  [0, 10)
    9  -- [0, 10)   [10, 20)
           ...       ...
    80 -- [80, 90)  [80, 90)
    """
    # TODO check input
    ra = to_deg(ra).value
    dec = to_deg(dec).value
    index_ra = np.digitize(ra, range(0, 100, 10)) - 1
    index_dec = np.digitize(dec, range(0, 100, 10)) - 1
    if (index_ra not in range(9)) or (index_dec not in range(9)):
        print("Error: ({}, {}) not in this octant.".format(ra, dec))
        patch_id = -1
    else:
        patch_id = index_dec * 9 + index_ra
    return patch_id


def which_patches(extent):
    """
    Determine which patches are covered by extent.

    Parameter
    ---------
    extent : float array, 4 values
        RA/Dec area given as [ramin, ramax, decmin, decmax].
    """
    # TODO check input
    ramin, ramax, decmin, decmax = extent
    p1 = which_patch(ramin, decmin) # lower left
    p2 = which_patch(ramax, decmin) # lower right
    p3 = which_patch(ramin, decmax) # upper left
    if not ((p1 >= 0) & (p2 >= 0) & (p3 >= 0)):
        patch_ids = []
    else:
        patch_ids = [range(y, y + p2 - p1 + 1) for y in range(p1, p3 + 9, 9)]
    return np.array(patch_ids).flatten()


def patch_extent(patch_id):
    """
    Determine RA/Dec boundaries of a Flagship patch.

    Parameter
    ---------
    patch_id : int in [0, 80]
        Flagship patch ID.
    """
    if patch_id not in range(81):
        print("Invalid patch_id.")
        extent = []
    else:
        ramin = (patch_id % 9) * 10.
        ramax = ramin + 10
        decmin = (patch_id / 9) * 10.
        decmax = decmin + 10
        extent = [ramin, ramax, decmin, decmax]
    return extent * u.deg


# def spv_extent(patch_id):
#     if patch_id not in (0, 1):
#         print("Invalid patch_id.")
#         return
#
#     if patch_id == 0:
#         extent = (220.5, 225.5, 12, 17)
#     else:
#         extent = (233, 238, 11.5, 16.5)
