# _Plasmodium vivax_ transmission model
This is the C++ source code for the _Plasmodium vivax_ transmission model developed by Dr. Michael WHITE (<mwhite@pasteur.fr>). This model has been first described in [White _et al._, 2018](dx.doi.org/10.1038/s41467-018-05860-8).

Briefly, it is a mixture of individual-based (humans) and compartmental models (mosquitoes), where _P. vivax_ infection in humans can lead to the introduction of hypnozoites (dormant stage of the parasite foudn in the liver). These hypnozoites then cause relapses through a stochastic process. The following figure summarises the different compartment.

![model](../resources/model_compartments.png)

This model also allows for simulating public health interventions. Currently implemented interventions are:
* LLINs: Distirbution of long-lasting insecticide nets, along with the progressive decay of the insecticide agent.
* IRS: Distribution and usage of insecticide-repellant spray, along with the progressive decay of insecticide agent.
* MDA (BS): Mass drug administration of bloodstage drug, e.g. artemisinin-combined therapy (ACT). This kind of MDA will clear bloodstage infections but will preserve the current hypnozoite reservoir.
* MDA (BS+PQ): Mass drug administration of both bloodstage (ACT) and hypnozoidicidal drug (Primaquine/Tafenoquine). This will clear ongoing bloodstage infections as well as dormant liver stages. Depending on options, administration of hypnozoidicidal drugs may be subject to restrictions (e.g. in G6PD-deficient individuals, pregnant women and infants).
* MSAT (BS+PQ): Mass screen-and-treat for current bloodstage infection, where people screened positive are treated with a bloodstage drug and a hypnozoidicidal drug, when it is safe to do so (_ie_ again exlucding G6PD-deficient individuals, pregnant women and infants)
* SSAT (BS+PQ): Serological screen-and-treat procedure, where the presence/absence of hypnozoites is assessed through the proxy of a recent bloodstage infection. This method uses ongoing developments of serological markers for a recent exposure to _P. vivax_. Individuals are screened and those who are thought the have had a bloodstage infections in the preceding 9 months are treated with hypnozoidicidal (again, when safety conditions are met). All individuals are given bloodstage drugs.

It should be noted that all these interventions can be tweaked with numerous parameters, such as  intervention coverage, drug efficiency, screening procedure (sensitivity and specifity).

## What's in this repo ?
There are two (main) branches to this GitLab repository :
* [Master](https://gitlab.pasteur.fr/mwhite/pv_mod/tree/master/), where the latest model developments are usually merged. This is the more up-to-date version of the model.
* [Legacy](https://gitlab.pasteur.fr/mwhite/pv_mod/tree/legacy), containing the version of the model that was initially developed and used in the abovementioned research paper. This branch is no longer actively maintained and much of the code has been refactored, especially in respect to treatment pathways.

The actual model code is found in the Pv\_mod folder, in the [Source.cpp](../Pv_mod/Source.cpp) file. Numerical optimisation routines are found in the same directory, in different C++ source and header files.
The root of this repository contains example of model parameter files, as well as some R code to visualize model outputs.

## Build the model
This section provides instructions on how to build the model on different operating systems.

### On your laptop
Under UNIX/Linux, use GCC 7+ (I know for sure build will fail with GCC 4 because of syntaxical errors with '\<' and '\>' signs).

```bash
g++-9 -O3 -o Pv_mod/Pv_model.o Pv_mod/Source.cpp Pv_mod/Simulation.cpp Pv_mod/Intervention.cpp Pv_mod/Population.cpp Pv_mod/Individual.cpp Pv_mod/Params.cpp Pv_mod/Mosquito.cpp Pv_mod/com.cpp Pv_mod/linpack.cpp Pv_mod/randlib.cpp 
```

TODO: write about building model under Windows

### On Pasteur TARS cluster (or any CentOS cluster)
You need to load GCC into your environment as it's not available by default.
```bash
module load gcc/8.2.0
g++ -O3 -o Pv_mod/Pv_model.o Pv_mod/Source.cpp Pv_mod/Simulation.cpp Pv_mod/Intervention.cpp Pv_mod/Population.cpp Pv_mod/Individual.cpp Pv_mod/Params.cpp Pv_mod/Mosquito.cpp Pv_mod/com.cpp Pv_mod/linpack.cpp Pv_mod/randlib.cpp 
```

## Use the model
Currently, arguments passed to the model binary are unnamed and should respect the following order:
1. General model parameters file
2. Peri-domestic mosquitos parameters file
3. Occupational mosqutos parameters file
5. Intervention parameters file
6. Model output file

Great care should be taken in generating these parameter files, especially the intervention file. 

### General model parameters
TODO: write about general model parameters

### Mosquito parameters
In former versions of the model, mosquito species were modeled explicitely and each had one input parameter file. We now model mosquito mosquitos regardless of their species, but assume they are found near houses (peri-domestic exposure) or at places were males are at increased exposure (occupational exposure). The fraction of males that are over-exposed is, for now, kept fixed at 37%, reflecting what can be observed in malaria-endemic countries. 
These two different mosquito populations make use of two input parameter files, [domestic_mosquitoes.txt](../domestic_mosquitoes.txt) and [occupational_mosquitoes.txt](../occupational_mosquitoes.txt).

### Intervention parameters
The [Intervention file generator](../Intervention_file_generator.R) R script will produce an intervention file where each column corresponds to a time where an intervention (or different combinations of interventions) is enforced.

It should be noted that this file has changed format between the [legacy branch](https://gitlab.pasteur.fr/mwhite/pv_mod/tree/legacy/Intervention_file_generator.R) and the [master branch](https://gitlab.pasteur.fr/mwhite/pv_mod/tree/master/Intervention_file_generator.R). The legacy version of the model used a row-based format, where each row corresponded to an intervention time. In the master branch, the newly-introduced format is column-based, as described in the previous paragraph. There are also more parameters, and not necessarily in the same order, in the master branch. 

__Therefore, intervention parameters files are not interchangeable.__

### Run the model
The following line will run the model using input files provided in this repository, assuming the model has been built in `Pv_mod/Pv_model.o` as described above:
```bash
Pv_mod/Pv_model.o model_parameters.txt domestic_mosquitoes.txt occupational_mosquitoes.txt intervention_coverage.txt model_output.txt
```
The output will be written to `model_output.txt`.

## Explore model output
Model output is provided in a text-based format, where each row represent a time step and columns corespond to model parameters. As such, the _i_-th line and _j_-th column will provide the value of parameter _j_ at time step _i_. The model simulates daily events, so that each line corresponds to a different day, with initial value of `365*date_start`.

The burn-in period (10 years) is dropped from model output, as well as some summary statistics that are not relevant for most uses. These can be further explored in  [Source.cpp, around line 195](../Pv_mod/Source.cpp#L195) (at the time of writing this ReadMe), where sections are commented out by default for the sake of making the output files a bit smaller.

### Summary statistics
TODO: write about `Output_overview.R` script