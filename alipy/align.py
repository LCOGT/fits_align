import star
import os
import numpy as np
import math
import scipy.ndimage
from astropy.io import fits
import csv

import logging

logger = logging.getLogger(__name__)



def affineremap(filepath, transform, shape, alifilepath=None, outdir = "alipy_out", hdu=0, verbose=True):
    """
    Apply the simple affine transform to the image and saves the result as FITS, without using pyraf.

    :param filepath: FITS file to align
    :type filepath: string

    :param transform: as returned e.g. by alipy.ident()
    :type transform: SimpleTransform object

    :param shape: Output shape (width, height)
    :type shape: tuple

    :param alifilepath: where to save the aligned image. If None, I will put it in the outdir directory.
    :type alifilepath: string

    :param hdu: The hdu of the fits file that you want me to use. 0 is primary. If multihdu, 1 is usually science.


    """
    inv = transform.inverse()
    (matrix, offset) = inv.matrixform()
    #print matrix, offset

    data, hdr = fromfits(filepath, hdu = hdu, verbose = verbose)
    data = scipy.ndimage.interpolation.affine_transform(data, matrix, offset=offset, output_shape = shape)

    basename = os.path.splitext(os.path.basename(filepath))[0]

    if alifilepath == None:
        alifilepath = os.path.join(outdir, basename + "_affineremap.fits")
    else:
        outdir = os.path.split(alifilepath)[0]
    if not os.path.isdir(outdir):
        os.makedirs(outdir)


    tofits(alifilepath, data, hdr = hdr, verbose = verbose)


def shape(filepath, hdu = 0, verbose=True):
    """
    Returns the 2D shape (width, height) of a FITS image.

    :param hdu: The hdu of the fits file that you want me to use. 0 is primary. If multihdu, 1 is usually science.


    """
    hdr = fits.getheader(filepath, hdu)
    if hdr["NAXIS"] != 2:
        raise RuntimeError("Hmm, this hdu is not a 2D image !")
    logger.debug("Image shape of %s : (%i, %i)" % (os.path.basename(filepath), int(hdr["NAXIS1"]), int(hdr["NAXIS2"])))
    return (int(hdr["NAXIS1"]), int(hdr["NAXIS2"]))


def fromfits(infilename, hdu = 0, verbose = True):
    """
    Reads a FITS file and returns a 2D numpy array of the data.
    Use hdu to specify which HDU you want (default = primary = 0)
    """

    logger.debug("Reading %s ..." % (os.path.basename(infilename)))

    hdul = fits.open(infilename)
    pixelarray = hdul[1].data
    # Why is this transpose done?
    pixelarray = np.asarray(pixelarray).transpose()
    hdr = hdul[1].header
    pixelarrayshape = pixelarray.shape
    logger.debug("FITS import (%i, %i) BITPIX %s / %s" % (pixelarrayshape[0], pixelarrayshape[1], hdr["BITPIX"], str(pixelarray.dtype.name)))

    return pixelarray, hdr

def tofits(outfilename, pixelarray, hdr = None, verbose = True):
    """
    Takes a 2D numpy array and write it into a FITS file.
    If you specify a header (fits format, as returned by fromfits()) it will be used for the image.
    You can give me boolean numpy arrays, I will convert them into 8 bit integers.
    """
    pixelarrayshape = pixelarray.shape
    logger.debug("FITS export (%i, %i) %s ..." % (pixelarrayshape[0], pixelarrayshape[1], str(pixelarray.dtype.name)))

    if pixelarray.dtype.name == "bool":
        pixelarray = np.cast["uint8"](pixelarray)

    if os.path.isfile(outfilename):
        os.remove(outfilename)

    if hdr == None: # then a minimal header will be created
        hdu = fits.PrimaryHDU(pixelarray.transpose())
    else: # this if else is probably not needed but anyway ...
        hdu = fits.PrimaryHDU(pixelarray.transpose(), hdr)

    hdu.writeto(outfilename)

    logger.debug("Wrote %s" % outfilename)
