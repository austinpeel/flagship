"""
This module contains methods for reading the Euclid Flagship Galaxy Mock.

_Author_ Austin Peel <austin.peel@cea.fr>

_Date_ 05 juillet 2017
"""

from common import home
import os
import astropy.io.fits as fits


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


def fetchcat(filename, *args):
    """
    Retrieve a Flagship galaxy catalog in fits format.

    Parameters
    ----------
    filename : str
        Name of the file as [***].fits

    Notes
    -----
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
