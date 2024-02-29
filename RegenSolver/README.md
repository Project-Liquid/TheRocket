## Regen Solver
This folder contains everything related to the solver developed for the Yale Project Liquid Team. Key files are:
1. <em>hemsolver.ipynb</em> - the main (HEM Model) solver
2. <em>solverscript.py</em> - a python script version of the solver

as well as some additional solver variations (with main differences related to coolant-property definition):

3. <em>subcomprsolver.ipynb</em> - a solver based on the main solver, but using properties for a coolant (using CoolProp) that is in the subcooled liquid state region
4. <em>subsupersolver.ipynb</em> - a solver based on the main solver, but using properties for a coolant (using Cantera) that is in the suberheated vapor state region

The main files to use, however are (1) and (2).
