"""
    FITS Align - Align and reproject FITS files from Las Cumbres Observatory
    Copyright (C) 2018 Edward Gomez

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import absolute_import
import sys, os
import math
import numpy as np
import operator
import copy
import itertools
from astropy.io import ascii

from fits_align.utils import cdist_np

import logging

logger = logging.getLogger(__name__)


class Star:
    """
    Simple class to represent a single source (usually stars, but not necessarily).
    In this module we often manipulate lists of such Star objects.
    """

    def __init__(self, x=0.0, y=0.0, name="untitled", flux=-1.0, props={}, fwhm=-1.0, elon=-1.0):
        """
        flux : Some "default" or "automatic" flux, might be a just good guess. Used for sorting etc.
        If you have several fluxes, colours, store them in the props dict.
        props : A placeholder dict to contain other properties of your choice (not required nor used by the methods).
        """
        self.x = float(x)
        self.y = float(y)
        self.name = str(name)
        self.flux = float(flux)
        self.props = props
        self.fwhm = float(fwhm)
        self.elon = float(elon)

    def copy(self):
        return copy.deepcopy(self)

    def __getitem__(self, key) :
        """
        Used for sorting list of stars.
        """
        if key == 'flux':
            return self.flux
        if key == 'fwhm':
            return self.fwhm
        if key == 'elon':
            return self.elon

    def __str__(self):
        """
        A string representation of a source.
        """
        return "%10s : (%8.2f,%8.2f) | %12.2f | %5.2f %5.2f" % (self.name, self.x, self.y, self.flux, self.fwhm, self.elon)

    def coords(self, full=False):
        """
        Returns the coords in form of an array.

        :param full: If True, I also include flux, fwhm, elon
        :type full: boolean

        """
        if full:
            return np.array([self.x, self.y, self.flux, self.fwhm, self.elon])
        else:
            return np.array([self.x, self.y])

    def distance(self, otherstar):
        """
        Returns the distance between the two stars.
        """
        return math.sqrt(np.sum((self.coords() - otherstar.coords())**2))

    def trigangle(self, otherstar):
        """
        returns the "trigonometric" angle of the vector that goes from
        self to the otherstar, in degrees
        """
        return math.atan2(otherstar.y - self.y, otherstar.x - self.x) * (180.0/math.pi) % 360.0

    def distanceandsort(self, otherstarlist):
        """
        Returns a list of dicts(star, dist, origpos), sorted by distance to self.
        The 0th star is the closest.

        otherstarlist is not modified.
        """
        import operator # for the sorting

        returnlist=[]
        for i, star in enumerate(otherstarlist):
            dist = self.distance(star)
            returnlist.append({'star':star, 'dist':dist, 'origpos':i})
        returnlist = sorted(returnlist, key=operator.itemgetter('dist')) # sort stars according to dist

        return returnlist

def listtoarray(starlist, full=False):
    """
    Transforms the starlist into a 2D numpy array for fast manipulations.
    First index is star, second index is x or y

    :param full: If True, I also include flux, fwhm, elon
    :type full: boolean

    """
    return np.array([star.coords(full=full) for star in starlist])


def area(starlist, border=0.01):
    """
    Returns the area covered by the stars.
    Border is relative to max-min
    """
    if len(starlist) == 0:
        return np.array([0, 1, 0, 1])

    if len(starlist) == 1:
        star = starlist[0]
        return np.array([star.x - 0.5, star.x + 0.5, star.y - 0.5, star.y + 0.5])

    a = listtoarray(starlist)
    (xmin, xmax) = (np.min(a[:,0]), np.max(a[:,0]))
    (ymin, ymax) = (np.min(a[:,1]), np.max(a[:,1]))
    xw = xmax - xmin
    yw = ymax - ymin
    xmin = xmin - border*xw
    xmax = xmax + border*xw
    ymin = ymin - border*yw
    ymax = ymax + border*yw
    return np.array([xmin, xmax, ymin, ymax])


def readmancat(mancatfilepath, verbose="True"):
    """
    Reads a "manual" star catalog -- by manual, I mean "not written by sextractor".
    So this is typically a *short* file.

    Comment lines start with #, blank lines are ignored.
    The format of a data line is

    starname xpos ypos [flux]

    The data is returned as a list of star objects.
    """

    if not os.path.isfile(mancatfilepath):
        logger.error("File does not exist :")
        logger.error(mancatfilepath)
        logger.error("Line format to write : starname xpos ypos [flux]")
        sys.exit(1)


    myfile = open(mancatfilepath, "r")
    lines = myfile.readlines()
    myfile.close

    table=[]
    knownnames = [] # We check for uniqueness of the names

    for i, line in enumerate(lines):
        if line[0] == '#' or len(line) < 4:
            continue
        elements = line.split()
        nbelements = len(elements)

        if nbelements != 3 and nbelements != 4:
            logger.error("Format error on line", i+1, "of :")
            logger.error(mancatfilepath)
            logger.error("The line looks like this :")
            logger.error(line)
            logger.error("... but we want : starname xpos ypos [flux]")
            sys.exit(1)

        name = elements[0]
        x = float(elements[1])
        y = float(elements[2])
        if nbelements == 4:
            flux = float(elements[3])
        else:
            flux = -1.0

        if name in knownnames:
            logger.error("Error in %s" % (mancatfilepath))
            logger.error("The name '%s' (line %i) is already taken." % (name, i+1))
            logger.error("This is insane, bye !")
            sys.exit(1)
        knownnames.append(name)

        #table.append({"name":name, "x":x, "y":y, "flux":flux})
        table.append(Star(x=x, y=y, name=name, flux=flux))


    logger.debug("I've read {} sources from {}".format(len(table), os.path.split(mancatfilepath)[1]))
    return table


def readsexcat(mycat):
    """
    :param cat: Catalogue from HDU[2]

    We read a sextractor catalog from hdu[2] and return a list of stars.
    Minimal fields that must be present in the catalog :

        * X
        * Y
        * FWHM
        * ELLIPTICITY
        * FLUX

    """
    returnlist = []

    # We check for the presence of required fields :
    minimalfields = ["x", "y", "fwhm", "ellipticity", "flux"]
    colnames = [x.lower() for x in mycat.names]
    for field in minimalfields:
        if field not in colnames:
            logger.error("Field %s not available in catalog extension of file !" % (field))
            sys.exit(1)

    logger.info("Number of sources in catalog : %i" % (len(mycat)))

    if len(mycat) == 0:
        logger.error("No stars in the catalog :-(")
    else :
        for i, mc in enumerate(mycat) :
            flux = mycat['FLUX'][i]
            fwhm = mycat['FWHM'][i]

            newstar = Star(x = mycat['X'][i], y = mycat['Y'][i], name = str(i), flux=flux,
                    fwhm = mycat['FWHM'][i], elon = mycat['ELLIPTICITY'][i])

            returnlist.append(newstar)

    logger.debug("I've selected %i sources" % (len(returnlist)))

    return returnlist

def findstar(starlist, nametofind):
    """
    Returns a list of stars for which name == nametofind
    """
    foundstars = []
    for source in starlist:
        if source.name == nametofind:
            foundstars.append(source)
    return foundstars

def sortstarlistbyflux(starlist):
    """
    We sort starlist according to flux : highest flux first !
    """
    sortedstarlist = sorted(starlist, key=operator.itemgetter('flux'))
    sortedstarlist.reverse()
    return sortedstarlist

def sortstarlistby(starlist, measure):
    """
    We sort starlist according to measure : lowest first !
    Where measure is one of flux, fwhm, elon
    """
    sortedstarlist = sorted(starlist, key=operator.itemgetter(measure))
    return sortedstarlist









class SimpleTransform:
    """
    Represents an affine transformation consisting of rotation, isotropic scaling, and shift.
    [x', y'] = [[a -b], [b a]] * [x, y] + [c d]
    """

    def __init__(self, v = (1, 0, 0, 0)):
        """
        v = (a, b, c, d)
        """
        self.v = np.asarray(v)

    def getscaling(self):
        return math.sqrt(self.v[0]*self.v[0] + self.v[1]*self.v[1])

    def getrotation(self):
        """
        The CCW rotation angle, in degrees
        """
        return math.atan2(self.v[1], self.v[0]) * (180.0/math.pi)# % 360.0

    def __str__(self):
        return "Rotation %+11.6f [deg], scale %8.6f" % (self.getrotation(), self.getscaling())


    def inverse(self):
        """
        Returns the inverse transform !
        """

        # To represent affine transformations with matrices, we can use homogeneous coordinates.
        homo = np.array([
        [self.v[0], -self.v[1], self.v[2]],
        [self.v[1],  self.v[0], self.v[3]],
        [0.0, 0.0, 1.0]
        ])

        inv = np.linalg.inv(homo)

        return SimpleTransform((inv[0,0], inv[1,0], inv[0,2], inv[1,2]))



    def matrixform(self):
        """
        Special output for scipy.ndimage.interpolation.affine_transform
        Returns (matrix, offset)
        """

        return (np.array([[self.v[0], -self.v[1]], [self.v[1], self.v[0]]]), self.v[2:4])


    def apply(self, x, y):
        """
        Applies the transform to a point (x, y)
        """
        xn = self.v[0]*x -self.v[1]*y + self.v[2]
        yn = self.v[1]*x +self.v[0]*y + self.v[3]
        return (xn, yn)

    def applystar(self, star):
        transstar = star.copy()
        (transstar.x, transstar.y) = self.apply(transstar.x, transstar.y)
        return transstar

    def applystarlist(self, starlist):
        return [self.applystar(star) for star in starlist]


def fitstars(uknstars, refstars, verbose=True):
    """
    I return the transform that puts the unknown stars (uknstars) onto the refstars.
    If you supply only two stars, this is using linalg.solve() -- perfect solution.
    If you supply more stars, we use linear least squares, i.e. minimize the 2D error.

    Formalism inspired by :
    http://math.stackexchange.com/questions/77462/
    """

    assert len(uknstars) == len(refstars)
    if len(uknstars) < 2:
        logger.debug("Sorry I cannot fit a transform on less than 2 stars.")
        return None

    # ukn * x = ref
    # x is the transform (a, b, c, d)

    ref = np.hstack(listtoarray(refstars)) # a 1D vector of lenth 2n

    uknlist = []
    for star in uknstars:
        uknlist.append([star.x, -star.y, 1, 0])
        uknlist.append([star.y, star.x, 0, 1])
    ukn = np.vstack(np.array(uknlist)) # a matrix

    if len(uknstars) == 2:
        trans = np.linalg.solve(ukn, ref)
    else:
        trans = np.linalg.lstsq(ukn, ref,rcond=None)[0]

    return SimpleTransform(np.asarray(trans))


def identify(uknstars, refstars, trans=None, r=5.0, getstars=False):
    """
    Allows to:
     * get the number or matches, i.e. evaluate the quality of the trans
     * get corresponding stars from both lists (without the transform applied)

    :param getstars: If True, I return two lists of corresponding stars, instead of just the number of matching stars
    :type getstars: boolean

    Inspired by the "formpairs" of alipy 1.0 ...

    """
    if trans != None:
        ukn = listtoarray(trans.applystarlist(uknstars))
    else:
        ukn = listtoarray(uknstars)
    ref = listtoarray(refstars)

    dists = cdist_np(ukn,ref) # Big table of distances between ukn and ref
    mindists = np.min(dists, axis=1) # For each ukn, the minimal distance
    minok = mindists <= r # booleans for each ukn
    minokindexes = np.argwhere(minok).flatten() # indexes of uknstars with matches

    logger.debug("%i/%i stars with distance < r = %.1f (mean %.1f, median %.1f, std %.1f)" % (np.sum(minok), len(uknstars), r,
            np.mean(mindists[minok]), np.median(mindists[minok]), np.std(mindists[minok])))

    matchuknstars = []
    matchrefstars = []

    for i in minokindexes: # we look for the second nearest ...
        sortedrefs = np.argsort(dists[i,:])
        firstdist = dists[i,sortedrefs[0]]
        seconddist = dists[i,sortedrefs[1]]
        if seconddist > 2.0*firstdist: # Then the situation is clear, we keep it.
            matchuknstars.append(uknstars[i])
            matchrefstars.append(refstars[sortedrefs[0]])
        else:
            pass # Then there is a companion, we skip it.

    del dists
    del ukn
    logger.debug("Filtered for companions, keeping %i/%i matches" % (len(matchuknstars), np.sum(minok)))

    if getstars==True:
        return (matchuknstars, matchrefstars)
    else:
        return len(matchuknstars)
