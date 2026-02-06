# Exploring the mechanisms and therapeutic strategies for anabolic resistance in sarcopenia: A systems modeling study

## Overview
This study builds upon the previously developed [McColl & Clarke (2024) Muscle Protein Synthesis Kinetic Model](https://www.cell.com/iscience/fulltext/S2589-0042(23)02711-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2589004223027116%3Fshowall%3Dtrue) to simulate potential causes of anabolic resistance observed in sarcopenia. 

We apply three complementary sensitivity analysis approaches to identify model parameters (i.e., physiological mechanisms) most likely to contribute to anabolic resistance:
1. **MPSA of the Kinetic Model** - A multi-parametric sensitivity analysis of the original McColl & Clarke (2024) model to identify model parameters most influencing muscle protein synthesis, muscle protein breakdown, and net balance.
2. **Naive Approach** - A broad, data-agnostic sensitivity analysis to identify parameters with the largest impact on anabolic resistance.
3. **Targeted Approach** - A data-informed simulation of putative mechanisms of anabolic resistance, as identified in the literature. Model parameters are varied based on literature-derived estimates to quantitatively assess their influence on muscle metabolism.

In the Targeted Approach, mechanisms are simulated both individually and in combination to evaluate how they may interact to reproduce the muscle metabolism responses obserbed in older adults.

This tool is being developed by the [Clarke Laboratory for Quantitative Exercise Biology](https://www.sfu.ca/clarkelab-bpk.html). Further information can be found in the following study.

> McColl TJ, Moore DR, Emberly E, Church DD, & Clarke DC. [Exploring the multifactorial causes and therapeutic strategies for anabolic resistance in sarcopenia: A systems modeling study.](https://www.biorxiv.org/content/10.1101/2025.09.12.675977v1) bioRxiv 2025.09.12.675977; doi: 10.1101/2025.09.12.675977

## Getting started
The 'McColl_2025_multifactorialCauseOfAnabolicResistance_250801' folder contains all required files to run the MPSA, simulations, and replicate the manuscript figures.

## Installation
Download the package to a local folder (e.g., '~/MultiFactorialCauseOfAnabolicResistance/') by extracting the ZIP file or by running the following terminal command:
```
git clone https://github.com/tjmccoll/MuscleProteinSynthesisKineticModel.git
```

### Running the simulations

The analyses presented in the paper are organized in three seperate folders:
1. KineticModelMPSA
2. NaiveApproach
3. TargetedApproach

Each folder contains a readMe file that details the included files in each folder

A seperate folder is included that contains the original McColl & Clarke (2024) _iScience_ model that underlies these analyses
* "McColl_2023_Muscle Protein Synthesis Kinetic Model_230919"

## Contact
tmccoll@sfu.ca

## Citation
If you use the simulations and analyses in this repository, please cite the paper:
> **McColl TJ, Moore DR, Emberly E, Church DD, & Clarke DC**\
> [Exploring the multifactorial causes and therapeutic strategies for anabolic resistance in sarcopenia: A systems modeling study.](https://www.biorxiv.org/content/10.1101/2025.09.12.675977v1)\
> bioRxiv 2025.09.12.675977\
> doi: 10.1101/2025.09.12.675977

If you use the mathematical model in this repository, please cite the original paper:
> **McColl, TJ & Clarke, DC.**\
> Kinetic modeling of leucine-mediated signaling and protein metabolism in human skeletal muscle\
> *iScience.* 27(1). 2024.\
> doi: [10.1016/j.isci.2023.108634](https://www.cell.com/iscience/fulltext/S2589-0042(23)02711-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2589004223027116%3Fshowall%3Dtrue)

