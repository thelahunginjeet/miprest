#!/usr/bin/env python

from distutils.core import setup,Command
from numpy.distutils.misc_util import get_numpy_include_dirs
from distutils.extension import Extension
import os

setup(name='miprest',
      version='1.1.0',
      description='Python Package for Mixed ICA/PCA',
      author='Kevin S Brown and Ameya Akkalkotkar',
      author_email='kevin.s.brown@uconn.edu',
      url='https://github.com/thelahunginjeet/miprest',
      packages=['miprest'],
      package_dir = {'miprest': ''},
      install_requires = ['pycar','pyica'],
      license='BSD-3',
      classifiers=[
          'License :: OSI Approved :: BSD-3 License',
          'Intended Audience :: Developers',
          'Intended Audience :: Scientists',
          'Programming Language :: Python',
          'Topic :: Blind Source Separation',
          'Topic :: Statistical Signal Processing',
      ],
    )
