# MSLibrarian
## Introduction 

In this workflow, we will use the R software package _MSLibrarian_ to create predicted spectral libraries for DIA proteomics analyses. MSLibrarian makes use of the DIA data at hand to incorporate calibrated predictions of both fragment ion intensities and retention times into the built spectral library. Apart from creating spectral libraries, it offers tools to manipulate existing spectral libraries on the protein, peptide and transition level. 

### Hardware requirements and recommedations

In its current form, MSLibrarian **must** be installed on a computer with **Windows** as the operating system. This requirement is mainly a consequence of the current third party softwares that the package uses for its operation. A future aim is to make MSLibrarian a cross-platform application as a docker image. 

A recommendation is to use a computer with at least 32 GB RAM to avoid issues during some of the more memory-requiring tasks that MSLibrarian performs. 

### Software requirements

To run MSLibrarian, the following softwares/pipelines must be installed on the **C:/**-drive: 

 * [**R version 4.0.0**](https://cran.r-project.org/) or higher. 
 * [**Proteowizard suite version 3.0.20365**](http://proteowizard.sourceforge.net/download.html) or higher 
 * [**Trans-proteomic pipeline version 5.2.0**](https://sourceforge.net/projects/sashimi/files/Trans-Proteomic%20Pipeline%20%28TPP%29/) or higher
 * [**OpenMS version 2.5.0**](https://github.com/OpenMS/OpenMS/releases/tag/Release2.6.0) or higher
 * [**PREGO**](https://bitbucket.org/searleb/prego-srm-response-predictor/downloads/) 
 * [**DeepLC GUI v0.1.29**](https://github.com/compomics/DeepLC/releases) or higher. _Follow the installation guide to setup the miniconda environment._

### Creating a predicted spectral library with _MSLibrarian_

### Download MSLibrarian from Github

To install _MSLibrarian_ from Github, you can make use of the **devtools** package. Currently, the repository is private, meaning that those with access needs to provide a [personal access token (PAT)](https://github.com/settings/tokens) to the argument _auth_token_ of the _install_github_ function. In the future, the repository will become a public one.  

```
library(devtools)
install_github("MarcIsak/MSLibrarian", auth_token = "change-to-your-personal-access-token")

```
