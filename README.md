# FITS Align

Align a sequence of astronomical FITS files based on sources extracted in each image. The aligned files will be geometrically reprojected so all the images are the same size and shape.

It is designed to work exclusively with reduced data from [Las Cumbres Observatory](https://lco.global). FITS data from LCO is [Rice compressed](https://heasarc.gsfc.nasa.gov/fitsio/fpack/) and contains sources extracted using [SEP](https://sep.readthedocs.io/en/v1.0.x/) during our data pipeline processing in a catalogue (`CAT`) FITS Header Data Unit (HDU).

## Installation

Download this repository :
`python setup.py install`

## Example usage

```python
from fits_align.ident import make_transforms
from fits_align.align import affineremap
from glob import glob
from numpy import shape

tmp_dir = "<FULL PATH TO INPUT FILES>"

img_list = sorted(glob(os.path.join(tmp_dir,"*.fz")))
ref_image = img_list[0]
images_to_align = img_list[1:]

identifications = make_transforms(ref_image, images_to_align)
```
If you have FITS files with image data you could use the following to reproject them.

```python
aligned_images = [ref_image]
for id in identifications:
    if id.ok:
        alignedimg = affineremap(id.ukn.filepath, id.trans, outdir=tmpdir)
        aligned_images.append(alignedimg)
```

If you just have FITS files which contain on photometry catalogues, which have pixel coordinate values (e.g. x, y) that you want to align (ie. no image data):

```python
import pandas as pd

identifications = make_transforms(img_list[0], img_list[1:], hdu='CAT')
catalogues = []
for i, catfile in enumerate(img_list):
    with fits.open(catfile) as hdul:
        data = pd.DataFrame.from_records(hdul['CAT'].data)
        if i != 0:
            (matrix, offset) = identifications[i-1].trans.matrixform()
            newcoords = np.dot(data[['x','y']], matrix) + offset
            (x,y) = np.transpose(newcoords)
            data.update({'x':x, 'y':y})
        catalogues.append(data)
```

## HDU containing catalogues

If you are using a none standard place for your catalogue HDU, that can be passed in a parameter to the `make_transforms()` function i.e.

```
identifications = make_transforms(ref_image, images_to_align, hdu=0)
```

## About

This is a customised fork of [AliPy by Malte Tewes](http://obswww.unige.ch/~tewes/alipy/). I have removed the dependency on SciPy in favour of pure NumPy (for linear algebra) and Pillow (for image array transforms).
