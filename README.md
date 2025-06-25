# ExBEMT: An Extended Blade Element Momentum Theory Solver

![image](https://github.com/user-attachments/assets/3ab02815-5df5-4a92-b018-612111ecc945)

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021b%2B-blue)](https://www.mathworks.com)

**ExBEMT (aka Extended-BEMT) is a modular BEMT-Based tool with some improvements for small-scale propellers and VTOL proprotors, optimized for high-incidence conditions and real-time simulation. Fast, accurate, and built for early-stage design and testing.**

---

## Introduction & Vision

**Version 1.0 (06/25/2025)**

Welcome to ExBEMT! My name is Ege Konuk, and this tool is a direct result of my PhD research in Aerospace Engineering at Old Dominion University. My work specializes in propeller performance prediction and developing reduced-order aerodynamic models for electric aircraft.

I believe that advanced engineering tools shouldn't be confined to commercial software or specialized research labs. My goal in making ExBEMT public is to provide a powerful, intuitive, and accessible open-source tool for students, hobbyists, and professional engineers. I hope it aids in your own research, design projects, and learning.

---

## Key Features

* **Interactive GUI:** A user-friendly graphical interface to easily set up and run complex analyses without modifying code.
* **Multiple Fidelity Models:** Seamlessly switch between different analysis models:
    * Standard Blade Element Momentum Theory (BEMT)
    * BEMT with a steady-state Pitt-Peters Dynamic Inflow model
    * BEMT with a Quasi-Unsteady Pitt-Peters Dynamic Inflow model
* **Automated Airfoil Analysis:** Integrates XFOIL to automatically generate or load airfoil polar data for different propeller sections.
* **Comprehensive Plotting:** Automatically generates performance curves (Thrust, Torque, Normal Force, Efficiency) and detailed azimuthal contour plots for in-depth analysis of blade loading.
* **Flexible Analysis Modes:** Analyze propeller performance over a range of advance ratios, velocities, or angles of incidence.
* **Modular Propeller Data:** Easily add and manage different propeller geometries through a simple folder structure.

---

## Showcase

### Main Application Interface

The intuitive GUI allows for complete control over the simulation parameters, from atmospheric conditions and RPM to the specific analysis range and fidelity model.

![Main ExBEMT GUI](https://raw.githubusercontent.com/your-username/your-repo/main/path/to/your/gui_screenshot.png)

### Detailed Azimuthal Contour Plots

Visualize the aerodynamic loads and flow conditions across the entire propeller disk for any flight condition. The tool automatically generates comparison plots to show the effects of different inflow models.

![Contour Plot Showcase](https://raw.githubusercontent.com/your-username/your-repo/main/path/to/your/contour_plot_screenshot.png)

---

## Analysis Models

ExBEMT implements several models to provide a range of fidelity and computational speed:

1.  **Baseline BEMT:** A robust and fast implementation of Blade Element Momentum Theory, suitable for initial performance estimates.
2.  **BEMT + Pitt-Peters Dynamic Inflow (Steady):** Incorporates the classic Pitt-Peters dynamic inflow model to capture the effects of non-uniform induced velocity, providing more accurate results for off-axis flow and transient conditions.
3.  **BEMT + Unsteady Pitt-Peters Dynamic Inflow:** The highest fidelity model, which implements the unsteady formulation of the Pitt-Peters model to capture time-varying induced velocity effects, making it ideal for high-incidence and complex VTOL flight regimes.

---

## Getting Started

### Prerequisites

* MATLAB (R2021b or newer)
* MATLAB Aerospace Toolbox (for `atmosisa` and `convangvel` functions)
* XFOIL: The executable (`xfoil.exe` for Windows) must be included in the project's root directory or in your system's PATH.

### Installation & Running

1.  **Clone the Repository:**
    ```bash
    git clone [https://github.com/your-username/ExBEMT.git](https://github.com/your-username/ExBEMT.git)
    ```
2.  **Open in MATLAB:** Navigate to the cloned repository folder in MATLAB.
3.  **Run the App:** Open the `ExBEMT_GUI.m` file in the editor and click the **Run** button, or type the following in the MATLAB command window:
    ```matlab
    ExBEMT_GUI
    ```

---

## How to Use

1.  **Select a Propeller:** Choose a propeller from the "Propeller" dropdown menu. (See "Adding New Propellers" below).
2.  **Set General Settings:** Input the flight altitude, number of blades, and RPM.
3.  **Choose Analysis Type:** Select whether to analyze over a range of Advance Ratio, Velocity, or Angle of Incidence.
4.  **Define Analysis Range:** Set the min, max, and step values for your chosen analysis type.
5.  **Select Fidelity Model:** Under "Corrections & Models", choose the desired analysis model.
6.  **Run Analysis:** Click the **"Run Analysis"** button to begin the simulation.
7.  **View Results:** Performance plots will appear in the "Output" panel. If "Plot Contour Plots" is checked, new figure windows will be generated for each iteration.

---

## Adding New Propellers

ExBEMT is designed to be modular. To add your own propeller, create a new folder inside the `Propeller Packages` directory. This folder must contain:

1.  **`Prop_sections.xlsx`:** An Excel file detailing the blade geometry. It should follow the format of the existing examples.
2.  **Airfoil `.dat` files:** Coordinate files for each airfoil used in the propeller design, with names corresponding to those in the `Prop_sections.xlsx` file.

Once the folder is created, restart the GUI, and the new propeller will appear in the dropdown menu.

---

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

---

## Citation

If you use ExBEMT in your research, please consider citing it:


Konuk, E. (2025). ExBEMT: An Extended Blade Element Momentum Theory Solver. GitHub. https://github.com/your-username/ExBEMT


Or, by citing the related research paper:

Konuk E., and Landman, D., "Computer Based Modeling for Tilt-Wing e-VTOL Propeller Performance," AIAA SciTech 2023 Forum, National Harbor, MD, 2023.


## Contact

For questions, issues, or collaboration, please open an issue on GitHub or contact Ege Konuk at [egekonuk@gmail.com](mailto:egekonuk@gmail.com).



![ExBEMT-Window](https://github.com/user-attachments/assets/3341133c-198a-443a-907d-3391ddf8106c)
