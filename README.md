 # Software and Processing Codes / Design Results

Welcome to the repository for the _Software & Processing Codes_ and _Design Results_ of the research project titled _**Viability Assessment of Large-Scale Claude Cycle Hydrogen Liquefaction: A Study on Technical and Economic Perspectives**_ by Panji B. Tamarona, Rene Pecnik, and Mahinder Ramdin.

## Contents

- `streamProps.py` contains a class called 'ThermoState' for obtaining the thermodynamic properties of streams from CoolProp (as a means to access REFPROP data).

### Folders

- **centrifugalCompressor:** The `compressor.py` file contains a class `Compressor` used for calculating centrifugal compressor parameters based on user input. The preliminary design procedure is implemented and executed in Jupyter Notebook. A walkthrough of the procedure is provided in `exampleComp.ipynb`. The compressor design results for the baseline and scale-up scenarios are also included.

- **radialInflowTurbine:** The `turbine.py` file contains a class `Turbine` used for calculating radial inflow turbine parameters based on user input. The preliminary design procedure is implemented and executed in Jupyter Notebook. A walkthrough of the procedure is provided in `exampleTurb.ipynb`. The turbine design results for the baseline and scale-up scenarios are also included.

- **plateFinHeatExchanger:** For non-catalytic Plate-Fin Heat Exchanger (PFHX) preliminary design, ASPEN EDR was utilized. The resulting EDR files are provided for reference. For catalytic PFHX, the preliminary design incorporates a combination of ASPEN EDR and a kinetic simulation tool developed by O'Neil et al. (2023). A walkthrough of the procedure is provided in `examplePFHX.ipynb`. To properly execute the code, it is necessary to install the hydrogen cryogenic PFHX chemical engineering model developed by O'Neil et al. (2023).

## Additional Resources

For more details and installation guide on the processing algorithm of the PFHX model by O'Neil et al. (2023), please refer to their original publication at: [https://doi.org/10.1016/j.cep.2023.109272](https://doi.org/10.1016/j.cep.2023.109272)

## Usage

If you download the entire repository as one folder, the code should work properly. We are currently working to make the code run in Google Colab.
