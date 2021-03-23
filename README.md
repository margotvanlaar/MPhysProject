# MPhysProject
Simulation code for MPhys Project: Modelling Private Medicine Outlets in Ghana


This simulation code models the distribution of pharmacies in Ghana, where the system is described by a statistical mechanics Canonical Enseble. A Monte Carlo Metropolis algorithm is used to simulate the evolution of the system.  

Pharmacies interact with their surroundings in two ways:
- They attract population within a radius, called the cut-off radius. This interaction is scaled with the population at that point.
- They repel other pharmacies in their proximity. The repulsion has a magnitude proportional to the area of overlap of two circles with a radius set by the attraction.


To run the code:

The main method is contained in the "MC_simulation.py" file. This imports three other modules: "dynamicsV4.py" which defines the dynamics of the system, "data_prep.py" which ensures data read into the code is in the right format, and "visualisation.py" which contains functions to visualise various aspects of the simulation.

The main method reads in a map of the population density distribution across Ghana, "gha_pd.xlsx".

The main method defines a number of system parameters at the start of the code, which can be altered by the used:
1. The number of pharmacies, n_pharma
2. The number of time steps performed by the Monte Carlo algorithm, n_sweeps
3. The attractive cut-off radius, r_c
4. The repulsive cut-off radius, r_r
5. The temperature of the system, T
6. The scale factor of the repulsive interaction, sf
