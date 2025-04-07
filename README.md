# RADAR-DESIGNER-MATLAB-TOOLBOX

## Overview
This repository contains MATLAB scripts and simulations for Radar Cross Section (RCS) modeling and analysis using MATLAB's Radar Designer. It includes scripts for loading RCS data, defining different RCS models, creating radar scenarios, and running simulations to evaluate target detection performance.

## Features
- Load and process RCS data from CSV files
- Define constant and frequency-dependent RCS models
- Set up and simulate a radar scenario in MATLAB
- Visualize RCS patterns and simulation results

## Repository Structure
- `Matlab_Rcs_Simulation.m`: Main script for setting up and running the radar simulation.
- `dummy_rcs_data.csv`: Example CSV file containing RCS values for testing.
- `results/`: Directory for storing simulation outputs and figures.

## Getting Started
### Prerequisites
- MATLAB (Tested with R2022b and later)
- Phased Array System Toolbox

### Usage
1. Clone the repository:
   ```bash
   git clone https://github.com/YOUR_USERNAME/RADAR-DESIGNER-MATLAB-TOOLBOX.git
   ```
2. Open MATLAB and navigate to the cloned directory.
3. Run `Matlab_Rcs_Simulation.m` to execute the simulation.
4. View the generated plots and results in the `results/` folder.

## Example Output
- Visualization of angle-dependent RCS
- Console logs showing simulation steps

## Future Work
- Implement more advanced RCS models
- Add support for different radar waveforms
- Enhance visualization with 3D plotting

## Authors

* **Velissarios Gkoulias** - *Initial work* - (https://github.com/VelissariosGkoulias)

## Extra Material
* https://www.mathworks.com/help/radar/ug/modeling-target-radar-cross-section.html
