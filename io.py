"""
This module contains methods for reading the Euclid Flagship Galaxy Mock
and halo catalogue.

Author : Austin Peel <austin.peel@cea.fr>

Date : 05 juillet 2017 (version 1.3.3)
       12 septembre 2017
       22 novembre 2017
       11 jan 2018 (updated to version 1.5.2)
"""
from common import home
import os
import astropy.io.fits as fits


def datapath(subdir=None):
    """
    Retrieve the system path to local Flagship files.

    It is assumed that files are located in $HOME/Data/Flagship/

    Parameters
    ----------
    subdir : str, optional
        Subdirectory(ies) to add after $HOME/Data/Flagship/

    Returns
    -------
    datapath : str
        Full system path to Flagship data files.
    """
    path = os.path.join(home, 'Data/Flagship')
    if subdir is not None:
        path = os.path.join(path, str(subdir))
    return path


def dboxpath(subdir=None):
    """
    Retrieve the system path to Dropbox Flagship files.

    It is assumed that files are located in $HOME/Dropbox/Data/Flagship/

    Parameters
    ----------
    subdir : str, optional
        Subdirectory(ies) to add after $HOME/Dropbox/Data/Flagship/
    """
    path = os.path.join(home, 'Dropbox/Data/Flagship')
    if subdir is not None:
        path = os.path.join(path, str(subdir))
    return path


def fetch_cat(patch_id, halos=False, version='1.5.2', verbose=False):
    """
    Retrieve a Flagship galaxy or halo catalog in fits format.

    Parameters
    ----------
    patch_id : int
        Catalogue ID number as .../patches/[patch_id]/galcat_full.fits
    halos : bool, optional
        Get halo catalog instead of galaxy catalog for this patch if True.
        Default is False.
    version : str, optional
        Flagship release version. Default is '1.5.2'.
    verbose : bool, optional
        Print loaded file path if True. Default is False.
    """
    if verbose:
        print("Flagship v" + version)

    if halos:
        halopath = '{}/patches/{}/halos/'.format(version, patch_id)
        # fpath = datapath(subdir=halopath, version=version)
        fpath = datapath(subdir=(halopath + 'halocat_13.5.fits'))
    else:
        catpath = '{}/patches/{}/'.format(version, patch_id)
        # fpath = datapath(subdir=catpath, version=version)
        fpath = datapath(subdir=(catpath + 'galcat_full.fits'))
    try:
        cat = fits.getdata(fpath)
        if verbose:
            print("Loaded " + fpath)
        return cat
    except IOError as e:
        print("Could not find {}".format(fpath))


def fetch_sandrines_patch():
    catpath = '1.3.3/patches/sandrine/Flagship_1736_cat_1.fits'
    return fits.getdata(datapath(catpath))


def fetch_spv(patch_id, halos=False, verbose=False):
    if halos:
        halopath = '1.3.3/patches/spv2/halos/halocat{}_13.5.fits'.format(patch_id)
        fpath = datapath(halopath)
    else:
        catpath = 'patches/spv2/spv2cat{}.fits'.format(patch_id)
        fpath = datapath(catpath, version='1.3.3')
    try:
        cat = fits.getdata(fpath)
        if verbose:
            print("Loaded " + fpath)
        return cat
    except IOError as e:
        print("Could not find {}".format(fpath))
        print("Might need to run split_spv.py first.")
