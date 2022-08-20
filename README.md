#NIU Undergrad Research: Granular Flow Simulation
LIGGGHTS + MatLab code for running granular flow simulations and processing/plotting relevent data.
## Installation steps

1.Please refer to GranularSimGuide.pdf for a full guide on required installation steps and maintenance

## Directory Walkthrough 

1. data : The data directory contains the pellet template currently being used by our multisphere simulation scripts
2. meshes : This directory contains all of the meshes for the conveyor belts and boundary wall.
3. PostProcessing : A directory containing MatLab scripts for determining pellet properties/dimensions, Sorting generated data, and Plotting that sorted data
4. scripts : A directory containing all of the LIGGGHTS simulation scripts that generate dump files required by the MatLab sorting scripts. in.fullbin is the main script.
