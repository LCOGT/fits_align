# About

*This is a fork of AliPy to GitHub to make my Docker install work. See link below to original*

This is a python package to quickly, automatically, and robustly identify geometrical transforms between optical astronomical images, using only field stars. The images can have different pixel sizes, orientations, pointings and filters.

http://obswww.unige.ch/~tewes/alipy/


# Installation

Quick :
`python setup.py install`

To create a source distribution :
`python setup.py sdist`

More info :
http://docs.python.org/distutils/


# To generate the documentation

```bash
cd doc
make apidoc
make html
```
Documentation is then available in:
`--> _build/index.html`
