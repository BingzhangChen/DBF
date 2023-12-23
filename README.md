# Disturbance and BEF README

This file accompanies the paper entitled "Environmental disturbance modulates the effect of diversity on productivity" by Chen submitted to Ecological Modelling and introduces each file (including the source data file) used for running models and statistical analysis as well as plotting.

## Authors

Bingzhang Chen

## Contact email

[bingzhang.chen\@strath.ac.uk](mailto:bingzhang.chen@strath.ac.uk){.email}

## Brief summary of the study

This study uses a simple phytoplankton-nutrient chemostat model to investigate how environmental disturbances (i.e., fluctuations in washout rate) affect the effect of phytoplankton diversity on productivity.

## Contributor

Bingzhang Chen is responsible for writing the Fortran and R code.

## LICENSE

All the codes are covered by the MIT license. Please see the **LICENSE** file for details.

# Metadata

## Software and packages

All the codes are written in R 4.2.0. The following R packages are used: **foreach** (version 1.5.2), **nlme** (version 3.1-157), **plyr** (version 1.8.7), **dplyr** (version 1.0.9). The  detailed package versions and environment can be found in the file **Session_Info.png**.  

## How to run the code

1.  Download the source code and data. The most updated code and data are available at <https://github.com/BingzhangChen/ActivationEnergy.git> which can be obtained either by using git clone or directly downloaded from Github.

2.  To run the main analysis (**Table 2** in the main paper), run **Table2.R** in R. This script also generates Fig. S1 and S2.

3.  To reproduce **Fig. 1** in the main paper, run **Fig1Eapp.R**.

4.  To reproduce **Fig. 2** in the main paper, run **Fig2.R**.

## Model source code

1. main.f90: The model source code written in Fortran 90. It needs to be compiled to generate an executable.

2. Makefile: the makefile used to compile the source code.

## Model output

1. 

## R scripts

1. NP_BEF.R: R script used for reading the model outputs and computing the average phytoplankton biomass and primary productivity at each environmental level. It also computes the diversity effects, selection, and complementarity effects for each environmental condition. 


## Other supplemental materials

1.  AmNat60700Suppl.pdf: Supplemental files including Supplement 1 (Derivations of Eq. 1 and Eq. 2 in the main text and additional analysis results of autotrophic and heterotrophic prokaryotes as well as insects (Table S1 and S2)) and Supplement 2 (Estimations of $E_{app}$ by incorporating cell size (Table S3)).





6. Session_Info.png: the output of the command **sessionInfo()** showing the detailed software package versions as well as dependencies.

## Funding

This study is supported by a Leverhulme Trust Research Project Grant (RPG-2020-389).