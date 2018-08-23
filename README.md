# FITS Align

Align a sequence of astronomical FITS files based on sources extracted in each image. The aligned files will be geometrically reprojected so all the images are the same size and shape.

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

identifications = make_transforms(ref_image, images_to_align)

aligned_images = [ref_image]
for id in identifications:
    if id.ok:
        alignedimg = affineremap(id.ukn.filepath, id.trans, outdir=tmpdir)
        aligned_images.append(alignedimg)

```

## About

This is a customised fork of [AliPy by Malte Tewes](http://obswww.unige.ch/~tewes/alipy/). I have removed the dependency on SciPy in favour of pure NumPy (for linear algebra) and Pillow (for image array transforms).
