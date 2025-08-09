# ExBEMT: An Extended Blade Element Momentum Theory Solver

![image](https://github.com/user-attachments/assets/3ab02815-5df5-4a92-b018-612111ecc945)

[![MATLAB](https://img.shields.io/badge/MATLAB-R2024a%2B-blue)](https://www.mathworks.com)

**ExBEMT (aka Extended-BEMT) is a modular BEMT-Based tool with some improvements for small-scale propellers and VTOL proprotors, optimized for high-incidence conditions and real-time simulation. Fast, accurate, and built for early-stage design and testing.**

---

## Introduction & Vision

**Version 3.0 (08/07/2025)**

Welcome to the next evolution of ExBEMT! My name is Ege Konuk, and this tool is a direct result of my PhD research in Aerospace Engineering at Old Dominion University. This major update marks a significant leap in fidelity and usability, directly integrating several key features from the future development roadmap of Version 2.0.

The headline feature is the integration of **NeuralFoil**, a deep-learning model for rapid and extensive aerodynamic coefficient prediction. This is complemented by a heavily upgraded **XFOIL workflow**, which now incorporates the **Viterna post-stall model**. This allows for the accurate extrapolation of airfoil data to high angles of attack (±180°), providing a complete and physically realistic drag polar, which is crucial for high-incidence VTOL flight regimes.

Furthermore, the data handling has been modernized. Polar data from XFOIL is now saved in **RPM-specific files** for higher accuracy, and both XFOIL and NeuralFoil results can be saved and loaded from efficient `.mat` files, streamlining the analysis process. The GUI has also been enhanced with live in-app plotting and a more robust data loading system that dynamically adapts to your propeller geometry.

I hope these powerful new features aid you in your own design projects, learning, and discovery!

---

## Key Features

* **Interactive GUI:** A user-friendly graphical interface to easily set up and run complex analyses without modifying code.
* **Advanced Airfoil Analysis Methods:**
    * **NeuralFoil Integration:** A new option to use a deep-learning model for generating comprehensive airfoil polars over a wide range of Reynolds numbers and angles of attack.
    * **Enhanced XFOIL with Viterna Post-Stall Model:** The XFOIL analysis pipeline now automatically blends results with the Viterna method to accurately model airfoil performance deep into the post-stall region.
* **Multiple Fidelity Models:** Seamlessly switch between different analysis models:
    * Standard Blade Element Momentum Theory (BEMT)
    * BEMT with a steady-state Pitt-Peters Dynamic Inflow model
    * BEMT with a Quasi-Unsteady Pitt-Peters Dynamic Inflow model
* **Flexible Data & Run Management:**
    * **RPM-Specific Polars:** XFOIL polar data is now generated and saved for specific RPMs, improving the accuracy of the simulation.
    * **`.mat` File Support:** Save and load airfoil polar data grids as `.mat` files for faster and more efficient analysis cycles.
* **Comprehensive Plotting:**
    * **Live In-App Performance Plots:** Automatically generates performance curves (Thrust, Torque, Normal Force, Efficiency) in real-time within the GUI tab for each run.
    * **Detailed Azimuthal Contour Plots:** Visualize aerodynamic loads and flow conditions across the entire propeller disk for any flight condition.
* **Flexible Analysis Modes:** Analyze propeller performance over a range of advance ratios, velocities, or angles of incidence.
* **Modular Propeller Data:** Easily add and manage different propeller geometries through a simple folder structure with dynamic section detection.

---

## Showcase

### Main Application Interface

The intuitive GUI allows for complete control over the simulation parameters, from atmospheric conditions and RPM to the specific analysis range, fidelity model, and polar generation method.

![ExBEMT_analysis_AOI](https://github.com/user-attachments/assets/5998a493-6020-41fc-b362-c273f3714392)

### Detailed Azimuthal Contour Plots

Visualize the aerodynamic loads and flow conditions across the entire propeller disk for any flight condition. The tool automatically generates comparison plots to show the effects of different inflow models.

![ExBEMT_Contours](https://github.com/user-attachments/assets/0156d7ad-16bb-47e4-99ee-b95519d3c28e)

---

## Analysis Models

ExBEMT implements several models to provide a range of fidelity and computational speed:

1.  **Baseline BEMT:** A robust and fast implementation of Blade Element Momentum Theory, suitable for initial performance estimates.
2.  **BEMT + Pitt-Peters Dynamic Inflow (Steady):** Incorporates the classic Pitt-Peters dynamic inflow model to capture the effects of non-uniform induced velocity, providing more accurate results for off-axis flow and transient conditions.
3.  **BEMT + Unsteady Pitt-Peters Dynamic Inflow:** The highest fidelity model, which implements the unsteady formulation of the Pitt-Peters model to capture time-varying induced velocity effects, making it ideal for high-incidence and complex VTOL flight regimes.

---

## Getting Started

### Prerequisites

* MATLAB (R2024a or newer)
* MATLAB Aerospace Toolbox (for `atmosisa` and `convangvel` functions)
* XFOIL: The executable (`xfoil.exe` for Windows) must be included in the project's root directory or in your system's PATH.
* **Python Environment (Optional):** Required for using the NeuralFoil analysis method. The solver will call the `run_neuralfoil.py` script.

### Installation & Running

1.  **Clone the Repository:**
    ```bash
    git clone [https://github.com/egekonuk/ExBEMT.git](https://github.com/egekonuk/ExBEMT.git)
    ```
2.  **Open in MATLAB:** Navigate to the cloned repository folder in MATLAB.
3.  **Run the App:** Open the `ExBEMT_GUI.m` file in the editor and click the **Run** button, or type the following in the MATLAB command window:
    ```matlab
    ExBEMT_GUI
    ```

---

## How to Use

1.  **Select a Propeller:** Choose a propeller from the "Propeller" dropdown menu.
2.  **Set General Settings:** Input the flight altitude, number of blades, RPM, and the desired azimuthal step size for the calculation.
3.  **Choose Analysis Type:** Select whether to analyze over a range of Advance Ratio, Velocity, or Angle of Incidence.
4.  **Define Analysis Range:** Set the min, max, and step values for your chosen analysis type.
5.  **Select Fidelity Model:** Under "Corrections & Fidelity", choose the desired BEMT analysis model.
6.  **Configure Run Options:**
    * Check **"Use Existing Polars"** to load previously generated data. This will activate the "Data Source" dropdown.
    * If generating new polars, select the **"Analysis Method"** (XFOIL or NeuralFoil).
    * Configure the settings for the chosen analysis method in the corresponding panel.
7.  **Run Analysis:** Click the **"Run Analysis"** button to begin the simulation.
8.  **View Results:** Live performance plots will appear in a new tab in the "Output" panel. If "Plot Contour Plots" is checked, new figure windows will be generated for each iteration.

---

## Adding New Propellers

ExBEMT is designed to be modular. To add your own propeller, create a new folder inside the `Propeller Packages` directory. This folder must contain:

1.  **`Prop_sections.xlsx`:** An Excel file detailing the blade geometry. The GUI now dynamically detects the number of sections, so you are not limited to a fixed number of columns.
2.  **Airfoil `.dat` files:** Coordinate files for each airfoil used in the propeller design, with names corresponding to those in the `Prop_sections.xlsx` file.
3.  **(Optional) Pre-computed Polar Data:** You can place RPM-specific `.xlsx` files or NeuralFoil `.mat` files in this folder. They will be automatically detected and made available in the "Data Source" dropdown when "Use Existing Polars" is checked.

Once the folder is created, restart the GUI, and the new propeller will appear in the dropdown menu.

---

## Acknowledgements

* The automated, parallel-processing interface with XFOIL was made possible thanks to the excellent `XFOILinterface` class developed by **Rafael Fernandes de Oliveira**. This tool is highly recommended for anyone looking to integrate XFOIL with MATLAB.
    * [theolivenbaum/XFOILinterface on GitHub](https://github.com/theolivenbaum/XFOILinterface)
* The deep-learning based airfoil analysis is powered by **NeuralFoil**, an excellent tool developed by **Peter Sharpe**.
    * [peterdsharpe/NeuralFoil on GitHub](https://github.com/peterdsharpe/NeuralFoil)

---

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](https://github.com/egekonuk/ExBEMT--An-Extended-Blade-Element-Momentum-Theory-Solver/blob/main/LICENSE) file for details.

---

## Citation

If you use ExBEMT in your research, please consider citing it:

Konuk, E. (2025). GitHub. https://github.com/egekonuk/ExBEMT

Or, by citing my PhD Dissertation:

Konuk E., "Computer Based Modeling for Small e-VTOL Propeller Performance," Ph.D Dissertation, Old Dominion University, 2024. https://digitalcommons.odu.edu/mae_etds/627/ DOI: 10.25777/8agh-jw47

## Contact

For questions, issues, or collaboration, please open an issue on GitHub or contact me :)
