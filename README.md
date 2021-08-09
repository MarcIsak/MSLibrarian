# MSLibrarian

## Requirements and recommendations

### Hardware 

In its current form, MSLibrarian **must** be installed on a computer with **Windows** as the operating system. This requirement is mainly a consequence of the current third party softwares that the package uses for its operation. A future aim is to make MSLibrarian a cross-platform application as a docker image. 

A recommendation is to use a computer with at least 32 GB RAM to avoid issues during some of the more memory-requiring tasks that MSLibrarian performs. 

### Software 

To run all features of MSLibrarian, the following softwares/pipelines **must** be installed on the **C:/**-drive

 * [**R version 4.0.0**](https://cran.r-project.org/) or higher. 
 * [**Proteowizard suite version 3.0.20365**](http://proteowizard.sourceforge.net/download.html) or higher 
 * [**Trans-proteomic pipeline version 5.2.0**](https://sourceforge.net/projects/sashimi/files/Trans-Proteomic%20Pipeline%20%28TPP%29/) or higher
 * [**OpenMS version 2.5.0**](https://github.com/OpenMS/OpenMS/releases/tag/Release2.6.0) or higher
 * [**PREGO**](https://bitbucket.org/searleb/prego-srm-response-predictor/downloads/) 
 * [**DeepLC GUI v0.1.29**](https://github.com/compomics/DeepLC/releases) or higher._Follow the installation guide to setup the miniconda environment._ 
  As an alternative, the **DeepLC CLI (.exe)** can be installed instead. 
 * [**DIA-NN v.1.8**](https://github.com/vdemichev/DiaNN/releases/tag/1.7.12). Currently the most recent version, but older version should also work. 

## Getting started 

### Download and install MSLibrarian from Github

To both download and install _MSLibrarian_ from Github, use the **devtools** package. 

```
library(devtools)
install_github("MarcIsak/MSLibrarian")

```
### Run MSLibrarian 

Go to the **WIKI** of this repository to learn how to create a predicted spectral library in MSLibrarian. There is also a page on how to manually create a Prosit prediction SQLite database. 


