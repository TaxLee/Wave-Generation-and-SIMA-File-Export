# Wave Generation and SIMA File Export Tool

**Author:** Shuijin LI  
**Email:** <lishuijin@nbu.edu.cn>  
**Date:** 2025.03.06

This MATLAB GUI tool allows users to generate wave time series and export the data as ASCII files in the SIMA format.

## Overview

The heart of this tool is the `waveGUI.m` file, which provides an intuitive interface for:

- Selecting the wave type (Regular or Irregular, with multiple spectrum models)
- Configuring simulation parameters such as time step, number of samples, ramp duration, and random seed
- Plotting the generated wave time series
- Exporting waves to SIMA-compatible ASCII files
- Saving and loading configuration settings

## Features

- **Wave Type Selection:** Toggle between Regular and Irregular waves.
- **Spectrum Models:** Choose from several options including JONSWAP, Torsethaugen, Ochi-Hubble, and PM variations.
- **Visualization:** Display the wave elevation in the GUI for real-time preview.
- **Configuration Management:** Quickly save or restore your simulation settings.
- **SIMA ASCII Export:** Easily create formatted text files for further analysis.

## Prerequisites

- MATLAB environment to run the GUI.
- Basic knowledge of MATLAB GUIs is beneficial.

## How to Use

1. **Launch the GUI:**
    - Run the `waveGUI.m` file in MATLAB by entering:

      ```
      >> waveGUI
      ```

2. **Set Simulation Parameters:**
    - Choose the wave type and corresponding spectrum model.
    - Define parameters such as time step, number of samples, ramp duration, etc.
3. **Preview the Wave:**
    - Click **Plot Wave** to visualize the generated wave.
4. **Export the Data:**
    - Click **Generate ASCII** to save the wave elevation in a SIMA-formatted ASCII file.
5. **Manage Configurations:**
    - Use **Save Config** and **Load Config** to handle your simulation settings.

## File Structure

- `waveGUI.m` – Main MATLAB file containing the GUI and related functions.
- `README.md` – This documentation file.
- Additional helper functions embedded within `waveGUI.m` for spectrum calculations, ramping, and file I/O operations.

## License

Please see the accompanying license file or documentation for details regarding usage rights and restrictions.
