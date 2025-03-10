# Wave Generation and SIMA File Export Tool

**Author:** Shuijin Li  
**Email:** <lishuijin@nbu.edu.cn>  
**Date:** 2025.03.06

This MATLAB GUI allows you to generate wave time series and export the data as SIMA-compatible ASCII files. It features an intuitive interface that caters to both beginners and advanced users, offering quick configuration and visualization of wave simulations.

## Overview

The core of the tool is the `waveGUI_v3.m` file which integrates robust spectrum modeling and seamless data export functionalities. The GUI supports:

- Toggle between Regular and Irregular wave types.
- Various spectrum models including multiple JONSWAP parameterizations, Torsethaugen (testing), Ochi-Hubble (testing), and PM variations.
- Customizable simulation parameters such as time step, number of samples, ramp duration, and random seed.
- Real-time plotting to preview generated wave elevations.
- Easy configuration management with options to load and save settings.

## Features

- **Wave Type & Spectrum Model:** Switch effortlessly between regular and irregular wave generation. The irregular wave mode provides a choice among several spectrum models, including enhanced formulations like Torsethaugen and Ochi-Hubble.
- **Simulation Parameters:** Define essential parameters (time step, sample count, ramp duration, etc.) to tailor the simulation to your specific requirements.
- **Visualization:** The GUI provides a built-in plotting area for a real-time preview of the wave elevation, enabling immediate feedback on parameter adjustments.
- **Data Export:** Export your generated wave data as a SIMA-formatted ASCII file. The updated file generation now produces a well-structured output:
   - **Line 1:** Number of data points (integer).
   - **Line 2:** Time step (floating point with 4 decimal places).
   - **Line 3:** A header comment line including tool information and author contact details.
   - **Line 4:** A single-line string listing the simulation parameters used.
   - **Subsequent Lines:** Wave elevation values, each formatted with a fixed `%12.6f` pattern.
- **Configuration Management:** Easily save your simulation settings and reload them later, ensuring quick setup changes and collaboration with colleagues.

## Prerequisites

- MATLAB with GUI support is necessary to run the `waveGUI_v3.m` file.
- Basic knowledge of MATLAB scripting and GUI interactions will help you get started.

## How to Use

1. **Launch the GUI:**
    - Run the `waveGUI_v3.m` file in MATLAB:  
       ``>> waveGUI_v3``
2. **Configure the Simulation:**
    - Select the desired wave type and corresponding spectrum model.
    - Adjust simulation parameters such as time step, number of samples, and ramp duration.
3. **Preview the Wave:**
    - Click **Plot Wave** to generate and view the wave elevation in real time.
4. **Export Data:**
    - Click **Generate ASCII** to export the wave data in SIMA-compatible ASCII format.
5. **Manage Configurations:**
    - Use the **Save Config** and **Load Config** options to store and retrieve your simulation settings easily.

## File Structure

- `waveGUI_v3.m` – Main MATLAB file containing the GUI and all related functions.
- `README.md` – This documentation file.
- Additional helper functions for spectrum calculations, data ramping, and file I/O operations are embedded within `waveGUI_v3.m`.

## License

This project is licensed under the MIT License. See the accompanying LICENSE file for detailed terms of use and restrictions.
