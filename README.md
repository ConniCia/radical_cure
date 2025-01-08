This is the source code for the following peer-reviewed publication

> Ciavarella C, Drakeley C, Price RN, Mueller I, White M. **Quantifying _Plasmodium vivax_ radical cure efficacy: a modelling study integrating clinical trial data and transmission dynamics**. The Lancet Infectious Diseases. 2025. https://doi.org/10.1016/S1473-3099(24)00689-3


## Coding credits

### _P. vivax_ recurrence model (PvRM) written in R

The R code for the PvRM has been developed by Constanze CIAVARELLA (constanze DOT ciavarella AT umontpellier DOT fr).


### _P. vivax_ individual-based model (PvIBM) written in C++

The C++ code for the PvIBM has been developed by Michael WHITE (michael DOT white AT pasteur DOT fr) anc collaborators, and is available as stand-alone code at https://gitlab.pasteur.fr/mwhite/pv_mod.

A copy of the C++ code (with some very minor edits that do NOT have an impact on mnodel dynamics) is included in this repository under [code/Pv_mod](code/Pv_mod).


### All other code

All other code has been developed by Constanze CIAVARELLA (constanze DOT ciavarella AT umontpellier DOT fr).


## Input data

This code uses input data from four sources:

- De-identified data on the time to first P vivax recurrence from the [IMPROV clinical trial](http://dx.doi.org/10.1016/S0140-6736(19)31285-1) is not publicly available and hence not included in this repository. Access to this data can be requested from the IMPROV study team. An empty placeholder file is included under [data/raw/placeholder_IMPROV_data.csv](data/raw/placeholder_IMPROV_data.csv) and read into the code in the script [code/R/input_A.R](code/R/input_A.R).
- Data from a meta-analysis of clinical trials by [Commons et al.](https://doi.org/10.1016/S1473-3099(23)00430-9) is included directly in the code in the script [code/R/input_B.R](code/R/input_B.R).
- Data from a meta-analysis of clinical trials by [Watson et al.](https://doi.org/10.7554/eLife.83433) is included directly in the code in the script [code/R/input_C.R](code/R/input_C.R).
- Posterior parameter draws obtained from the calibration of the PvIBM by [White et al](http://dx.doi.org/10.1038/s41467-018-05860-8) are included in input file [data/raw/PvMod_posterior_draws.csv](data/raw/PvMod_posterior_draws.csv).


## Getting started

This code should be run **interactively** from the [code/main.R](code/main.R) script. It is thus recommended to use RStudio.

We use a cluster to run computationally heavy parts of the code. All cluster operations are conveniently handled through the R interface using ssh connections. You can specify the cluster address and the path on the cluster in the [code/R/_header.R](code/R/_header.R) file.

R packages to install are 'ssh' and those listed in [code/R/_header.R](code/R/_header.R).

C++ code is compiled using gcc version 9.2.0.


## Repository structure

- [code/main.R](code/main.R) is the R interface from where to run analysis steps.
- [code/R/](code/R/) is the directory containing R scripts.
- [code/Pv_mod/](code/Pv_mod/) is the directory containing the C++ code.
- [code/bash/](code/bash/) is the directory containing the bash files used to run cluster jobs.
- [data/raw/](data/raw/) is the directory containing raw input data.

Additional directories and files will be created within the [data/](data/) directory when running the code.


## Code structure

The code is comprised of four **analysis** steps (**A**, **B**, **C**, **D**) and a section to create visualisations of analysis results (diagnostic MCMC plots, figures and tables).

The code should be run **in sequence** since later parts of the analysis use output from earlier parts.


### Estimating primaquine and tafenoquine efficacy from clinical trial data

The PvRM simulates recurrent _P. vivax_ blood-stage infections in trial participants treated with primaquine or tafenoquine for a symptomatic _P. vivax_ infection. To estimate the efficacy of a primaquine or tafenoquine regimen, we calibrate the PvRM to clinical trial data via Markov Chain Monte Carlo (MCMC).

In code **step A**, we fit the PvRM to data from the [IMPROV clinical trial](http://dx.doi.org/10.1016/S0140-6736(19)31285-1) to estimate the hypnozoiticidal efficacy of 7 mg/kg of primaquine over 7 days and of 7 mg/kg of primaquine over 14 days.

In code **step B**, we fit the PvRM to data from a meta-analysis of clinical trials by [Commons et al.](https://doi.org/10.1016/S1473-3099(23)00430-9) to estimate the hypnozoiticidal efficacy of 3.5 mg/kg of primaquine over 7 or 14 days.

In code **step C**, we fit the PvRM to data from a meta-analysis of clinical trials by [Watson et al.](https://doi.org/10.7554/eLife.83433) to estimate the hypnozoiticidal efficacy of 5 mg/kg of tafenoquine and of 7.5 mg/kg of tafenoquine.


### Quantifying the community-level impact of case management with primaquine or tafenoquine

The PvIBM describes the dynamics of _P. vivax_ transmission at the population level. In particular, the PvIBM simulates _P. vivax_ relapses, G6PD activity levels by sex, case management, public health interventions, population demography, heterogeneity in exposure to mosquito bites, antiparasite and antidisease immunity, mosquito seasonality, larval mosquito stages, and vector control.

In code **step D**, we use the PvIBM to simulate the introduction of primaquine or tafenoquine in routine case management of symptomatic _P. vivax_ infections.


## Abbreviations used in the code

General abbreviations

-   8AQ = 8-Aminoquinoline (class of antibiotic drugs comprehending PQ and TQ)
-   BS = blood stage
-   CQ = chloroquine (BS drug)
-   DP = dihydroartemisinin-piperaquine (BS drug)
-   MCMC = Markov Chain Monte Carlo
-   LS = liver stage
-   PQ = primaquine (LS drug)
-   Pv = _Plasmodium vivax_ (parasite causing recurring malaria)
-   TQ = tafenoquine (LS drug)


Drug regimen abbreviations

-   BS_CQ: CQ over 3 days
-   BS_DP: DP (exact treatment course not specified)
-   PQ_lowdose7 :  \~ 3.5 mg/kg PQ total dose over 7 days (\~Brazil)
-   PQ_highdose7 : \~ 7 mg/kg PQ total dose over 7 days
-   PQ_lowdose14 : \~ 3.5 mg/kg PQ total dose over 14 days
-   PQ_highdose14: \~ 7 mg/kg PQ total dose over 14 days
-   TQ_lowdose : \~ 5 mg/kg TQ single dose (300 mg TQ for 60 kg person)
-   TQ_highdose: \~ 7.5 mg/kg TQ single dose (450 mg TQ for 60 kg person)
-   AQ_ideal: a hypothetical single-dose 8AQ drug
