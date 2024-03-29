# Fassa
Geometric VOF interface advection using an Eulerian Implicit method (Weymouth and Yue) and field operations on cartesian grids.

## Val di Fassa
Fassa does not stand for Fast Adaptive Simulator Something... (it is not adaptive neither fast). Instead, it is named after *Val di Fassa*, which is a beautiful valley located in the middle of Italian Dolomites. It is also where I wrote this code, which was simply intended as a school project.

## Geometric VOF Interface Advection
In CFD, the Volume Of Fluid (VOF) approach can be used to describe a two-phase system, by defining a scalar marker function ``f`` which represents the volume fraction of a reference phase in each computational cell of the domain. This volume fraction takes value of 1 in the pure phase, 0 in the other phase, and intermediate values at the interface between the two phases.
Advecting the interface according to the velocity field can be achieved from the solution of the following advection equation:
```math
\frac{\partial f}{\partial t} + u\cdot\nabla(f)=0
```  
Fassa solves this equation using a geometrical discretization technique of the advection term, in order to avoid interface smearing. The implemented approach was proposed by Weymouth and Yue, it requires a uniform cartesian grid, and it is conservative if a divergence-free velocity field ``u`` is provided.
> Weymouth, Gabriel D., and Dick K-P. Yue. "Conservative volume-of-fluid method for free-surface simulations on cartesian-grids." Journal of Computational Physics 229.8 (2010): 2853-2865

## Requirements
Fassa is written in python3 and it just requires numpy.
The volume fraction fields are initialized using the [Vofi library](https://github.com/VOFTracking/Vofi) which must be installed only if different configurations needs to be simulated.
The simulation results are .vtk files that can be visualized using [Paraview](https://www.paraview.org).

## Installation using 'pip'
```sh
$ git clone https://github.com/edocipriano/fassa.git
$ cd fassa
$ python3 setup.py bdist_wheel
$ pip3 install dist/fassa-1.1.0-cp38-cp38-macosx_11_0_x86_64.whl
```

The path of the last instruction may change based on the OS, and the version of the package.

## Structure
Fassa is organized in three main folders: paris, fassa, run.
* **fassa**: the main module of this project. It contains a small collection of classes and functions that allow to easily peform VOF simulations, to solve the pressure-velocity coupling using a fairly simple Projection Method, and to write the simulation results in vtk format.
* **paris**: a module that contains a few functions taken from the [Paris simulator](http://www.ida.upmc.fr/~zaleski/paris/) code. These functions are written in fortran and compiled using f2py in order to be imported from the python code. This module performs geometric operations, such as finding the interface normal using the Young's method, finding the plane constant, and cutting a cubic cell with a plane in order to find the fraction of cell occupied by the reference phase.
> Aniszewski, Wojciech, et al. "Parallel, robust, interface simulator (PARIS)." Computer Physics Communications 263 (2021): 107849.
* **run**: contains examples of simulations performed using fassa. The most remarkable cases are **vortex** and **zalesak**, which are classical VOF-advection benchmark simulations. In these simulations, a circle is advected by a prescribed velocity field and the final configuration is compared with the initial position of the circle, in order to test the accuracy of the implemented advection method.
> Tryggvason, Grétar, Ruben Scardovelli, and Stéphane Zaleski. Direct numerical simulations of gas–liquid multiphase flows. Cambridge university press, 2011.

### Run
After the installation process, one can simply run the examples typing:
```sh
$ cd run/vortex
$ python3 vortex.py
```
The simulation will create a folder called **results** where the vtk results are stored. After some post-processing using Paraview, the vortex results will look like this:

<img src="./run/vortex/fig-vortex.png" alt="vortex" width="1200"/>
