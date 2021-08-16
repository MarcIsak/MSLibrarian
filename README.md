# MSLibrarian

## Requirements and recommendations

### Hardware 

In its current form, MSLibrarian **must** be installed on a computer with **Windows** as the operating system. This requirement is mainly a consequence of the current third party softwares that the package uses for its operation. A future aim is to make MSLibrarian into a cross-platform application, and also provide it as a docker image. 

A recommendation is to use a computer with at least 32 GB RAM to avoid issues during some of the more memory-requiring tasks that MSLibrarian performs. 

### Software 

To run all features of MSLibrarian, the following softwares/pipelines **must** be installed on the **C:/**-drive

 * [**R version 4.0.0**](https://cran.r-project.org/) or later. 
 * [**Proteowizard suite version 3.0.20365**](http://proteowizard.sourceforge.net/download.html) or later 
 * [**Trans-proteomic pipeline version 5.2.0**](https://sourceforge.net/projects/sashimi/files/Trans-Proteomic%20Pipeline%20%28TPP%29/) or later
 * [**MSFragger version 3.2**](http://msfragger-upgrader.nesvilab.org/upgrader/) or later 
 * [**OpenMS version 2.5.0**](https://github.com/OpenMS/OpenMS/releases/tag/Release2.6.0) or later
 * [**PREGO**](https://bitbucket.org/searleb/prego-srm-response-predictor/downloads/) 
 * [**DeepLC GUI version 0.1.29**](https://github.com/compomics/DeepLC/releases) or later._Follow the installation guide to setup the miniconda environment._ 
  As an alternative, the **DeepLC CLI (.exe)** can be installed instead. 
 * [**DIA-NN version 1.8**](https://github.com/vdemichev/DiaNN/releases/tag/1.7.12). Currently the most recent version, but older versions should also work. 

## MS Data 

The input MS data must conform to the following: 

* **Format:** Thermo raw (other file formats should be available in the future) 
* **Acquisition mode:** DIA (must contain both MS1 and MS2 scans)
 
## Getting started 

### Download and install MSLibrarian from Github

To both download and install _MSLibrarian_ from Github, use the [**devtools**](https://cran.r-project.org/web/packages/devtools/index.html) package. 

```
library(devtools)
install_github("MarcIsak/MSLibrarian")

```
### Prosit prediction SQLite databases

MSLibrarian relies on Prosit prediction SQLite databases to make spectral libraries. SQLite databases can be downloaded for some 
common species from ProteomeXchange(**insert PXD submission**), or be manually created in MSLibrarian. Go to the **Wiki** of this repository to learn how 
to create your own Prosit prediction SQLite database. 

### Parameter files 

Some of the third-party tools that MSLibrarian uses, such as Comet (or MSFragger), DIA-Umpire or Spectrast require parameters files. Example parameter files can be found in the folder **params** in this repository. Make sure to add the folder to the C:/-drive and that the path to the files **do not** contain any spaces. 

These parameter files can be edited, but it is **not recommended** to do so.

### Run MSLibrarian 

Go to the **Wiki** of this repository to learn how to create a predicted spectral library in MSLibrarian. 


