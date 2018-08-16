# LCO-AliPy

This is a python package to quickly, automatically, and robustly identify geometrical transforms between optical astronomical images, using only field stars. The images can have different pixel sizes, orientations, pointings and filters.

It is designed to work exclusively with reduced data from [Las Cumbres Observatory](https://lco.global). FITS data from LCO is [Rice compressed](https://heasarc.gsfc.nasa.gov/fitsio/fpack/) and contains sources extracted using [SEP](https://sep.readthedocs.io/en/v1.0.x/) during our data pipeline processing.

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
hdu = fits.open(ref_image)
data = hdu[1].data
outputshape = shape(data)

identifications = make_transforms(ref_image, images_to_align)

for id in identifications:
    if id.ok:
        affineremap(id.ukn.filepath, id.trans, shape=(outputshape[1],outputshape[0]), outdir=tmpdir)

aligned_images = sorted(glob(tmpdir+"/*_affineremap.fits"))
```

## About

This is a customised fork of [AliPy by Malte Tewes](http://obswww.unige.ch/~tewes/alipy/). I have removed the dependency on SciPy in favour of pure NumPy (for linear algebra) and Pillow (for image array transforms).
