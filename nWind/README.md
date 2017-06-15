# **Timothy Stovall Neutrally Stable ABL - Reproduction Attempt**

----------
This repository contains the OpenFOAM case files used in my first attempt to reproduce the results of Timothy Stovall's neutrally stable atmospheric boundary layer (ABL) within his master's thesis.

----------


Overview
-------------
There are some differences between the basic setup of my simulation and that of Timothy Stovall.  Mr. Stovall was persuing a simulation which accurately represented the energy deficit within the wake of wind turbines, using actuator disk theory.  My persuit was to accurately represent a neutrally stable ABL at lower average wind speeds using the methods of Mr. Stovall.  It is also worth mentioning that OpenFOAM-3.0.0 was used for my simulation.

### Computational Domain
The computational domain in both simulations was the same, where L<sub>x</sub> x L<sub>y</sub> x L<sub>z</sub> = 1250 m x 800 m x 800 m in the streamwise, vertical, and spanwise directions, respectively.  The horizontal directions, in both meshes, had a mesh spacing of 10 m.  However, our vertical mesh grading method was a little different.  I chose to use an overall expansion ratio of 10:1 from the ground to the top boundary.  Mr. Stovall, however, started with a 2.5 m spacing at the ground.  At 40 m the cell height increased at a rate of 1.2, until a cell height of 10 m at an elevation of 90 m.  From 90 m to 800 m, the cell height remained a constant 10 m.

### Initial Conditions
Mr. Stovall initialized his vertical velocity profile using a power law, where *alpha*=0.097, *y_0*=100 m, and *U_0*=10 m/s.  My simulation was not initialized with a velocity profile.

### Boundary Conditions
Periodic boundary conditions were applied in the streamwise direction, while symmetry conditions were applied at the spanwise boudaries.  For velocity, a symmetry condition was applied at the top and a no-slip condition was applied at the ground.  For pressure, a symmetry condition was applied at the top and a zero-gradient condition was applied at the ground.  This matched Mr. Stovall's boundary conditions.

### Solver
Mr. Stovall utilized the *oodles* solver, which has been replaced by the *pimpleFoam* and is what I used.  Stovall used the one-equation eddy viscosity LES model (*OneEqEddy*).  I used various LES models, with and without the use of van Driest damping.


### Momentum Source
A momentum source is used to ensure the flow contains enough energy to develop and maintain turbulence.  My simulation utilizes the *fvOption* called *meanVelocityForce*, which adjusts the pressure field in order to maintain a specified bulk velocity (3 m/s, in this case).  To my knowledge, stovall did not use a momentum source.  However, he did use "roughness blocks" to increase turbulence intensity.

References
----------
- [Stovall, T., 2009, "Simulations of Wind Turbine Wake Interactions in OpenFOAM," M.S. thesis, Dept. Mech. Eng., Univ. CO.][1]

[1]: http://gradworks.umi.com/14/69/1469025.html