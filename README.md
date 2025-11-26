![](images/logo.png)

![](https://cranlogs.r-pkg.org/badges/grand-total/Coxmos)

* [Introduction](https://github.com/BiostatOmics/Coxmos#introduction)

* [Installation](https://github.com/BiostatOmics/Coxmos#installation)

* [Getting started](https://github.com/BiostatOmics/Coxmos#getting-started)

* [Contact](https://github.com/BiostatOmics/Coxmos#contact)

* [References](https://github.com/BiostatOmics/Coxmos#references)


### Introduction
The **Coxmos** R package is an end-to-end pipeline designed for the study of survival analysis for 
high dimensional data. Updating classical methods and adding new ones based on sPLS technologies. 
Furthermore, includes multiblock functions to work with multiple sets of information to improve 
survival accuracy. 

The pipeline includes three basic analysis blocks:

1. Computing cross-validation functions and getting the models. 

2. Evaluating all the models to select the better one for multiple metrics.

3. Understanding the results in terms of the global model and the original variables.

*Coxmos* contains the necessary functions and documentation to obtain from raw data the final models
after compare them, evaluate with test data, study the performance individually and in terms of 
components and graph all the results to understand which variables are more relevant for each case 
of study.

### Installation

#### Dependencies requiring manual installation

Some of the metrics available in *Coxmos* are optional based and will not be included in the 
standard *Coxmos* installation. A list of all optional packages are shown below:

* nsROC:
* smoothROCtime:
* survivalROC:
* risksetROC:
* ggforce:
* RColorConesa:

#### Installing Coxmos

The *Coxmos* R package and all the remaining dependencies can be installed 
from CRAN:

```
install.packages("Coxmos")
```

Or from GitHub using `devtools`

```
devtools::install_github("BiostatOmics/Coxmos")
```

In case of using Github, to access vignettes, you will need to force building with
`devtools::install_github(build_vignettes = TRUE)`. Please note that this will
also install all suggested packages required for vignette build and might 
increase install time. Alternatively, an HTML version of the vignette is
available under the [vignettes](https://github.com/BiostatOmics/Coxmos/tree/master/vignettes)
folder.


### Getting started

In order to use Coxmos, you will need the following items:

- An explanatory X matrix.
- A response survival Y matrix (with two columns, "time" and "event").

Please note that two toy datasets are included in the package. Details to load and use them can be 
found in the package's vignette.


### Contact
If you encounter a problem, please 
[open an issue](https://github.com/BiostatOmics/Coxmos/issues) via GitHub.

  
### References
If you use *Coxmos* in your research, please cite the original publication:

Coxmos: Interpretable survival models for high-dimensional and multi-omic data
Pedro Salguero, Anabel Buendía-Galera, Sonia Tarazona
bioRxiv 2025.10.01.679727; doi: [https://doi.org/10.1101/2025.10.01.679727](https://doi.org/10.1101/2025.10.01.679727)
