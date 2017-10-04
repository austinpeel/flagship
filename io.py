"""
This module contains methods for reading the Euclid Flagship Galaxy Mock
and halo catalogue.

Author : Austin Peel <austin.peel@cea.fr>

Date : 05 juillet 2017
       12 septembre 2017
"""
from common import home
import os
import astropy.io.fits as fits
# from .utils import make_kappa


def datapath(*args):
    """
    Retrieve the system path to the catalogs.

    It is assumed that the catalogs are fits files located in
    $HOME/Data/Flagship/

    Parameters
    ----------
    args : str or str list
        Name(s) of the subdirectory(ies) to add after $HOME/Data/Flagship/
    """

    path = os.path.join(home, 'Data/Flagship')
    for arg in args:
        path = os.path.join(path, arg)

    return path


def fetch_cat(filename, *args):
    """
    Retrieve a Flagship galaxy or halo catalog in fits format.

    Parameters
    ----------
    filename : str
        Name of the file as [***].fits

    Note
    ----
    Directories specifying the path to `filename` can be given as arguments
    after `filename`. For example,
    ```python
    fetchcat('cat.fits', 'subdir', 'split5', 'tiles')
    ```
    is equivalent to
    ```python
    fetchcat('subdir/split5/tiles/cat.fits')
    ```
    """

    filepath = os.path.join(datapath(*args), filename)
    if os.path.exists(filepath):
        cat = fits.getdata(filepath)
        return cat
    else:
        print("Could not find {}".format(filepath))


def fetch_kappa(split_id=0, zbin=None, npix=1024):
    """
    Retrieve a Flagship kappa map.

    Parameters
    ----------
    split_id : int
    zbin : int
    npix : int
    """
    if zbin is None:
        zb = ''
    else:
        zb = zbin
    filename = 'splits/{}/maps/kappa{}_{}.fits'.format(split_id, zb, npix)
    filepath = datapath(filename)
    try:
        kappa = fits.getdata(filepath)
    except IOError:
        start = len(datapath())
        print("...{} does not exist. Generating.".format(filepath[start:]))
        # if zbin is None:
        #     make_kappa(split_id=split_id, npix=npix)
        #     return fetch_kappa(split_id=split_id, npix=npix)
        # else:
        #     make_kappa(split_id, zbin, npix)
        #     return fetch_kappa(split_id, zbin, npix)
        print("Just kidding. Generate yourself.")

    return kappa
