# Wave Generation and SIMA File Export Tool

**Author:** Shuijin Li  
**Email:** <lishuijin@nbu.edu.cn>  
**Date:** 2025.03.06

This MATLAB GUI enables users to generate wave time series while seamlessly exporting the data to SIMA-compatible ASCII files. Its intuitive interface allows both beginners and advanced users to configure and visualize wave simulations quickly.

## Overview

At the core of the tool is the `waveGUI.m` file, which combines robust spectrum modeling and data export functionalities. The GUI supports:

- Selection between Regular and Irregular wave types.
- Multiple spectrum models including various JONSWAP parameters and PM variations.
- Custom simulation parameters: time step, sample count, ramp duration, and random seed.
- Real-time preview of generated wave elevation.
- Configuration management with option to save and load settings.

## Features

- **Wave Type & Spectrum Model:**  
    Easily switch between Regular and Irregular wave generation. For irregular waves, choose from models such as JONSWAP (different parameter sets) and PM variations.

- **Simulation Parameters:**  
    Configure essential values like time step, number of samples, and ramp duration to tailor the simulation to your needs.

- **Visualization:**  
    A built-in plotting area provides a real-time preview of the wave elevation, ensuring quick feedback on parameter adjustments.

- **Data Export:**  
    Export the generated wave data as a SIMA-formatted ASCII file, perfect for further analysis or integration into other tools.

- **Configuration Management:**  
    Use the Load/Save configuration options to quickly restore previous simulation settings or share configurations with colleagues.

## Prerequisites

- MATLAB (with GUI support) is required to run the `waveGUI.m` file.
- Familiarity with MATLAB scripting and basic GUI interactions is helpful.

## How to Use

1. **Launch the GUI:**
     - Run the `waveGUI.m` file in MATLAB:
         ```
         >> waveGUI
         ```

2. **Configure the Simulation:**
     - Select the wave type and corresponding spectrum model.
     - Set parameters for the simulation including time step, number of samples, ramp duration, etc.

3. **Preview the Wave:**
     - Click **Plot Wave** to generate and visualize the wave elevation in the GUI.

4. **Export Data:**
     - Use the **Generate ASCII** button to produce a SIMA-formatted ASCII file containing the wave data.

5. **Manage Configurations:**
     - Save your settings using **Save Config** and reload them later with **Load Config** for quick setup changes.

## File Structure

- `waveGUI.m` – Main MATLAB file including the GUI and related functions.
- `README.md` – This updated documentation file.
- Additional helper functions for spectrum calculations, data ramping, and file I/O operations are embedded within `waveGUI.m`.

## License

This project is licensed under the MIT License. See the accompanying LICENSE file for detailed usage rights and restrictions.
