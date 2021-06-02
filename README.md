# MSLibrarian

## Requirements and recommendations

### Hardware 

In its current form, MSLibrarian **must** be installed on a computer with **Windows** as the operating system. This requirement is mainly a consequence of the current third party softwares that the package uses for its operation. A future aim is to make MSLibrarian a cross-platform application as a docker image. 

A recommendation is to use a computer with at least 32 GB RAM to avoid issues during some of the more memory-requiring tasks that MSLibrarian performs. 

### Software 

<span style="color:blue">To run MSLibrarian, the following softwares/pipelines **must** be installed on the **C:/**-drive</span>: 

 * [**R version 4.0.0**](https://cran.r-project.org/) or higher. 
 * [**Proteowizard suite version 3.0.20365**](http://proteowizard.sourceforge.net/download.html) or higher 
 * [**Trans-proteomic pipeline version 5.2.0**](https://sourceforge.net/projects/sashimi/files/Trans-Proteomic%20Pipeline%20%28TPP%29/) or higher
 * [**OpenMS version 2.5.0**](https://github.com/OpenMS/OpenMS/releases/tag/Release2.6.0) or higher
 * [**PREGO**](https://bitbucket.org/searleb/prego-srm-response-predictor/downloads/) 
 * [**DeepLC GUI v0.1.29**](https://github.com/compomics/DeepLC/releases) or higher. _Follow the installation guide to setup the miniconda environment._
 * [**DIA-NN v.17.12**](https://github.com/vdemichev/DiaNN/releases/tag/1.7.12). Currently the most recent version, but older version should also work. 

## Getting started 

### Download and install MSLibrarian from Github

To both download and install _MSLibrarian_ from Github, use the **devtools** package. Currently, the repository is private, meaning that those with access needs to provide a [personal access token (PAT)](https://github.com/settings/tokens) to the argument _auth_token_ of the _install_github_ function. In the future, the repository will become a public.  

```
library(devtools)
install_github("MarcIsak/MSLibrarian", auth_token = "change-to-your-personal-access-token")

```
### Run MSLibrarian 

To use MSLibrarian from your library folder, you simply run: 

```
library(MSLibrarian)

```
To learn how to run MSLibrarian, go the subfolder **examples** of this repository and follow the tutorial described in **demo_mslibr.rmd**.


