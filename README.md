# LCO-AliPy

This is a python package to quickly, automatically, and robustly identify geometrical transforms between optical astronomical images, using only field stars. The images can have different pixel sizes, orientations, pointings and filters.

It is designed to work exclusively with reduced data from [Las Cumbres Observatory](https://lco.global). FITS data from LCO is [Rice compressed](https://heasarc.gsfc.nasa.gov/fitsio/fpack/) and contains sources extracted using [SEP](https://sep.readthedocs.io/en/v1.0.x/) during our data pipeline processing.

## Installation

Download this repository :
`python setup.py install`


## To generate the documentation

```bash
cd doc
make apidoc
make html
```
Documentation is then available in:
`--> _build/index.html`

## About

This is a customised fork of [AliPy by Malte Tewes](http://obswww.unige.ch/~tewes/alipy/). 
