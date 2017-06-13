# **NIST Technical Note 1944 - Reproduction Attempt #1**

----------
This repository contains the OpenFOAM case files used in my first attempt to reproduce the results of Dr. Liang Shi and Dr. DongHun Yeo within the NIST Technical Note 1944.

----------


Overview
-------------
There are many differences between the basic setup of my simulation and that of Technical Note 1944.  With the study having been done using custom tools and solvers, I attempted to recreate their results using standard OpenFOAM utilities and solvers.  It is also worth mentioning that OpenFOAM-3.0.0 was used for my research, while Dr. Shi and Dr. Yeo used OpenFOAM-2.3.0.

### Computational Domain
The computational domain is at 1:100 scale of the atmospheric boundary layer (ABL), where L<sub>x</sub> x L<sub>y</sub> x L<sub>z</sub> = 20 m x 10 m x 20 m in the streamwise, vertical, and spanwise directions, respectively.  The grid contains cells that are uniform and equidistant in every direction, resulting in 150 x 75 x 150 cells.  Dr. Shi and Dr. Yeo set up their domain with the y-direction as the spanwise direction and the z-direction as the vertical direction.  I chose to swap these to be more consistent with other simulations within my research.

### Initial Conditions
Dr. Shi and Dr. Yeo created a utility, based off the work of Schoppa and Hussain, for "kickstarting" turbulence.  I attempted to utilize a similar utility created by Eugene de Villers, called perturbU, however the mesh resolution near the ground was much too course for it to have any effect.  As a result, I opted to let turbulence develop on its own, but initialized the flow with a parabolic vertical velocity profile.

### Boundary Conditions
Periodic boundary conditions are applied in all horizontal directions for all quantities.  For velocity, a slip condition is applied at the top and a no-slip condition is applied at the ground.  For pressure, a zero-gradient condition is applied at both the top and the ground.  These conditions do not vary from the NIST paper.  However, the Dr. Shi and Dr. Yeo utilized a wall function from Schumann.  There is no wall function applied in my simulation.

### Solver
The solver created for the NIST paper is based off *pimpleFOAM*, a solver packaged with OpenFOAM.  However, extra features were added by Dr. Shi and Dr. Yeo for additional statistical processing.  They also used the one-equation eddy viscosity SGS model, *kEqn/oneEqEddy*, but implemented a custom wall function based off Schumann LES model.  I started the simulation with no wall function, but implemented the *nutkWallFunction* after having difficulty developing turbulence.

### Momentum Source
A momentum source is used to ensure the flow contains enough energy to develop and maintain turbulence.  Dr. Shi and Dr. Yeo created a custom *fvOption* called *constGradPForce*, which maintains a constant pressure gradient along the streamwise direction.  My simulation utilizes the *fvOption* called *meanVelocityForce*, which adjusts the pressure field in order to maintain a specified bulk velocity (27 m/s, in this case).

The *fvSolution* file was adjusted to match the NIST paper specifications, where time discretization is performed by the implicit Crank-Nicolson scheme, using a coefficient of 0.5.  Pressure is solved by the generalized geometric-algebraic multi-grid (GAMG) algorithm.  Velocity is solved by the preconditioned bi-conjugate gradient solver for asymmetric matrices (PBiCG) algorithm, using the simplified diagonal-based incomplete LU (DILU) preconditioner.  The use of DILU was not explicitly mentioned in the NIST paper.
