# Deviated Fixed-route Microtransit
Tools to model microtransit network design and operations routing, used for the modeling and computations of the paper "Deviated Fixed-route Microtransit: Design and Operations".

# Usage
This repository contains code written using Julia v1.5.


- The Main.jl file contains information about the required Julia packages and the content of each of the code files.
- The data directory contains an input folder containing the data to run the multiple algorithms, and an output folder containing the trial information of the computational experiments conducted in the paper.
- The pipeline files start with "MiND" and there is one for each direct model of the VRP variant (path, subpath, and segment-based), one for each Benders-based algorithmic method (variable enumeration, and exact and heuristic column generation), and the subpath-based implementation of the DAR variant.
- To run a trial, you can do so using the command window and running the following command from the directory of the pipeline file: 
  
  `julia "pipeline_file" -o "path to the output directory" -f "name of the trial info file" -t "trial id"`.
- It is recommended to place the trial info file in the same location as the desired output directory.