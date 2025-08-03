# Exploring the mechanisms and therapeutic strategies for anabolic resistance in sarcopenia: a computational approach

## Overview
This study builds upon the previously developed [McColl & Clarke (2024) Muscle Protein Synthesis Kinetic Model](https://www.cell.com/iscience/fulltext/S2589-0042(23)02711-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2589004223027116%3Fshowall%3Dtrue) to simulate potential causes of anabolic resistance observed in sarcopenia. 

We apply two complementary sensitivity analysis approaches to identify model parameters (i.e., physiological mechanisms) most likely to contribute to anabolic resistance:
1. **Naive Approach** - A broad, data-agnostic sensitivity analysis to identify parameters with the largest impact on anabolic resistance.
2. **Targeted Approach** - A data-informed simulation of putative mechanisms of anabolic resistance, as identified in the literature. Model parameters are varied based on literature-derived estimates to quantitatively assess their influence on muscle metabolism.

In the Targeted Approach, mechanisms are simulated both individually and in combination to evaluate how they may interact to reproduce the muscle metabolism responses obserbed in older adults.

This tool is being developed by the [Clarke Laboratory for Quantitative Exercise Biology](https://www.sfu.ca/clarkelab-bpk.html). Further information can be found in the following study.

> (update with BioRxIv link)

## Getting started
The 'McColl_2025_multifactorialCauseOfAnabolicResistance_250801' folder contains all required files () to run the MPSA, simulations, and replicate the manuscript figures.

## Installation
Download the package to a local folder (e.g., '~/MultiFactorialCauseOfAnabolicResistance/') by extracting the ZIP file or by running the following terminal command:
```
git clone https://github.com/tjmccoll/MuscleProteinSynthesisKineticModel.git
```

### Running the simulations

To simulate the study analyses:
1. ...
2. ...

### File list 

**(UPDATE)**

#### 'Code' folder:
* ‘McColl_2023_model execution_230922.mlx’: This script contains the code to simulate the model and create all the plots included in the manuscript and supplementary information. 
#### 'Experimental Data' folder:
* ‘230221_experimental data.xlsx’: This excel file contains all experimental data collected for model calibration and validation.
* ‘200832_Biolo, 3-pool parameters.xlsx’: 3-pool parameter values at baseline. These experimental values were used to calibrate the baseline k-values that control the 3-pool parameters. The comparison between the simulated and experimental baseline 3-pool parameter values are calculated at lines 42-43 in the ‘McColl_2023_model execution_230922.mlx’ script.
* '230308_Unit conversions.xslx': An excel file that provides the calculations for the conversion of units from experimental data sets to the input data in the model.
#### 'Functions' folder:
* Contains all required functions to run the model execution script. Information pertaining to each function is contained within the function script.
#### 'Input Data' folder:
* ’230922_Initial Values, IOM.xlsx’: Contains the initial values for each species that is inputted into the ODE function. 
* ’231023_K-values, IOM.xlsx’: Contains the kinetic rate parameters that are inputted into the ODE function.
#### 'Meta-analysis' folder:
* Contains the R code and experimental data needed replicate the p-Akt and p-p70S6K meta-analysis and spline regressions.
#### 'Optimizer' folder:
* Contains the ‘runOptimizer_230923’ code to run the parameter optimizer and the ‘objectiveFunction_230923’ function that the optimizer seeks to minimize.
#### 'Output Plots' folder:
* Empty folder
* The model execution script will create a folder corresponding to the date and time that line 9 is run where saved model figures are stored.

## Contact
tmccoll@sfu.ca

## Citation
If you use the simulations and analyses in this repository, please cite the paper:
> **(ADD CITATION INFO)**

If you use the mathematical model in this repository, please cite the original paper:
> **McColl, TJ & Clarke, DC.**\
> Kinetic modeling of leucine-mediated signaling and protein metabolism in human skeletal muscle\
> *iScience.* 27(1). 2024.\
> doi: [10.1016/j.isci.2023.108634](https://www.cell.com/iscience/fulltext/S2589-0042(23)02711-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2589004223027116%3Fshowall%3Dtrue)

