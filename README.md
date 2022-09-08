# heater_membrane
## Introduction
Here is an example simulation, showing the steady-state thermal response of a 2-D system. A wire heater is pattered on a silicon-dioxide membrane, as might be manufactured with a silicon-on-insulator process. For example, this may be used in a flow sensor. An advection term is used in the PDE so that the flow of water, or air, can be considered.

## Directory structure
-src<br/>
pdeSim.py - code for solving the PDE<br/>
models.py - contains the model geometry<br/>
materials.csv - thermal properties of materials, referenced by models.py<br/>

-scripts<br/>
membrane_heater.ipynb - runs the overall simulation, performs data analysis<br/>

-data<br/>
membrane_heater.hdf5 - the simulation data, including the input parameters, as attributes<br/>

-figures<br/>
*.png - snapshots of the simulation results<br/>

## Methodology
The PDE is written as a sparse matrix equation and directly solved using the scipy.sparse.linalg.spsolve method. At a larger matrix size, an iterative solved would have been used but, here, the run time of the direct solver is acceptable.