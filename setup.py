#!/usr/bin/env python
from setuptools import setup

setup(
    name='fits_align',
    version='0.4.4',
    packages=['fits_align'],
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    description="Align a sequence of astronomical FITS files based on sources extracted in each image. The aligned files will be geometrically reprojected so all the images are the same size and shape.",
    author='Edward Gomez',
    install_requires=[
        'astropy',
        'numpy',
        'pillow'
    ],

    classifiers=[
    "Development Status :: 5 - Production/Stable",
    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    'Topic :: Software Development :: Build Tools',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
],
)
