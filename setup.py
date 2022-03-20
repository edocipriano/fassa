import os
from setuptools import find_packages
from numpy.distutils.core import setup, Extension

description = ("fassa: Geometric VOF interface advection using an Eulerian Implicit method (Weymouth and Yue) and field operations on cartesian grids.\n")

ext1 = Extension(name='paris.al3d',
                sources=['paris/al3d.pyf',
                         'paris/al3d.f90'])

ext2 = Extension(name='paris.fl3d',
                sources=['paris/fl3d.pyf',
                         'paris/fl3d.f90'])

ext3 = Extension(name='paris.area3d',
                sources=['paris/area3d.pyf',
                         'paris/area3d.f90'])

setup(name='fassa',
      packages=find_packages(),
      package_dir={'fassa': 'fassa'},
      version='1.1.0',
      description='fassa',
      keywords='navier-stokes; fvm; vof; advector',
      url='https://github.com/edocipriano/fassa',
      author='Edoardo Cipriano',
      author_email='edoardo.cipriano@polimi.it',
      license='MIT',
      install_requires=['numpy','scipy','numba'],
      test_suite='',
      ext_modules = [ext1, ext2, ext3]
    )
