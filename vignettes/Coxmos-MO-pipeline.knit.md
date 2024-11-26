---
title: "Step-by-step guide to the MO-Coxmos pipeline"
author: 
    name: "Pedro Salguero García"
    affiliation: "Polytechnic University of Valencia and Institute for Integrative Systems Biology (I2SysBio), Valencia, Spain"
    email: pedrosalguerog[at]gmail.com
package: Coxmos
date: "Last updated: 11 noviembre, 2024"
output: 
    rmarkdown::html_vignette:
        toc: true
vignette: >
  %\VignetteIndexEntry{Step-by-step guide to the MO-Coxmos pipeline}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

<style type="text/css">
.main-container {
  max-width: 1080px !important;
  margin-left: auto;
  margin-right: auto;
}

body {
  max-width: 1080px !important;
  margin-left: auto;
  margin-right: auto;
}
</style>



# Introduction

The **Coxmos** R package includes the following basic analysis blocks:

1. **Cross-validation and Modeling for High-Dimensional Data**: The user can use the *Coxmos* package to select the optimal parameters for survival models in high-dimensional datasets. The package provides tools for estimating the best values for the parameters.

2. **Comparing Classical and High-Dimensional Survival Models**: After obtaining multiple survival models, the user can compare them to determine which one gives the best results. The package includes several functions for comparing the models.

3. **Interpreting Results**: After selecting the best model or models, the user can interpret the results. The package includes several functions for understanding the impact of the original variables on survival prediction, even when working with (s)PLS methods.

4. **Predicting New Patients**: Finally, if a new dataset of patients is available, the user can use the model to make predictions for the new patients and compare the variables against the model coefficients to estimate the patients' risk of an event.

# Installation

Coxmos can be installed from GitHub using `devtools`:


```r
install.packages("devtools")
devtools::install_github("BiostatOmics/Coxmos", build_vignettes = TRUE)
```

# Getting ready

To run the analyses in this vignette, you'll first need to load `Coxmos`:


```r
# load Coxmos
library(Coxmos)
```

In addition, we'll require some additional packages for data formatting. Most of them are signaled as `Coxmos` dependencies, 
so they will already be installed in your system.

To generate plots, we make use of the `RColorConesa` R package.
After install:


```r
# install.packages("RColorConesa")
library(RColorConesa)
```

# Input data

The Coxmos pipeline requires two matrices as input. First one must be the features under study, and the second one a matrix with two columns called **time** and **event** for survival information.

After quality control, the data contains expression data for miRNAs and Protein information. The data could be load as follows:


```r
# load dataset
data("X_multiomic", package = "Coxmos")
data("Y_multiomic", package = "Coxmos")

X <- X_multiomic
Y <- Y_multiomic

rm(X_multiomic, Y_multiomic)
```

These files contain a `list` of two blocks and a `data.frame` object. Our toy example has a total of **150 observations, 642 miRNAs and 369 proteins**:



|                 | hsa-let-7a-2-3p| hsa-let-7a-3p| hsa-let-7a-5p| hsa-let-7b-3p| hsa-let-7b-5p|
|:----------------|---------------:|-------------:|-------------:|-------------:|-------------:|
|TCGA-A2-A0SV-01A |        9.344533|     228.00660|      56661.51|      155.1192|      57564.19|
|TCGA-A2-A0YT-01A |       25.031847|     134.54618|      56991.26|      103.2564|      64478.91|
|TCGA-BH-A1F0-01A |       34.725664|     122.98673|      68110.05|      195.3319|      37871.23|
|TCGA-B6-A0IK-01A |       44.962897|     163.66495|      89585.87|      183.4486|     121455.58|
|TCGA-E2-A1LE-01A |       15.392269|      87.54353|     188649.57|      109.6699|      54892.68|



|                 |     7529|      7531|      7534|      1978|     7158|
|:----------------|--------:|---------:|---------:|---------:|--------:|
|TCGA-A2-A0SV-01A | 0.025273| -0.030661|  0.458350|  1.369000| -0.40781|
|TCGA-A2-A0YT-01A | 0.348210|  0.347870|  0.891570| -0.188540| -0.13896|
|TCGA-BH-A1F0-01A | 0.132340| -0.348210|  0.078923|  0.227890| -0.49063|
|TCGA-B6-A0IK-01A | 0.352670|  0.376670| -0.285750| -0.380980| -1.43080|
|TCGA-E2-A1LE-01A | 0.110140|  0.094990| -0.159960| -0.098148|  0.40160|



|                 | time|event |
|:----------------|----:|:-----|
|TCGA-A2-A0SV-01A |  825|TRUE  |
|TCGA-A2-A0YT-01A |  723|TRUE  |
|TCGA-BH-A1F0-01A |  785|TRUE  |
|TCGA-B6-A0IK-01A |  571|TRUE  |
|TCGA-E2-A1LE-01A |  879|TRUE  |



As can be observed, clinical variables were transform to binary/dummy variables for factors.


```r
ggp_density.event <- plot_events(Y = Y, 
                                 categories = c("Censored","Death"), #name for FALSE/0 (Censored) and TRUE/1 (Event)
                                 y.text = "Number of observations", 
                                 roundTo = 0.5, 
                                 max.breaks = 15) 
```


```r
ggp_density.event$plot
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-147-1.png" width="100%" />

# Survival models for multiomic data sets

After loading the data, it may be of interest for the user to perform a survival analysis in order to examine the relationship between explanatory variables and the outcome. However, traditional methods are only applicable for low-dimensional datasets. To address this issue, we have developed a set of functions that make use of (s)PLS techniques in combination with Cox analysis for the analysis of high-dimensional datasets.

*Coxmos* provides the following methodologies for multi-omic approaches:

* SB.sPLS approaches: SB.sPLS-ICOX, iSB.sPLS-ICOX, SB.sPLS-DRCOX-Penalty, iSB.sPLS-DRCOX-Penalty, SB.sPLS-DRCOX-Dynamic, iSB.sPLS-DRCOX-Dynamic, SB.sPLS-DACOX-Dynamic and iSB.sPLS-DACOX-Dynamic.
* MB.sPLS approaches: MB.sPLS-DRCOX-Dynamic and MB.sPLS-DACOX-Dynamic.

More information for each approach could be found in the help section for each function. The function name for each methodology are:

* SB.sPLS approaches: `sb.splsicox()`, `isb.splsicox()`, `sb.splsdrcox()`, `isb.splsdrcox()`, `sb.splsdrcox()`, `isb.splsdrcox()`, `sb.splsdacox()` and `isb.splsdacox()`.
* MB.sPLS approaches: `mb.splsdrcox()` and `mb.splsdacox()`.

To perform a survival analysis with our example, we will use methodologies that can work with high-dimensional data. These are the set of methodologies that use sPLS techniques.

The first thing we are going to do is split our data into a train and test set. This split will be made with a proportion of 70% of the data for training and 30% for testing.

## TRAIN/TEST

We will use the function `createDataPartition` from the R package **caret**. We will use a 70% - 30% split for training and testing, respectively, and set a seed for reproducible results.


```r
set.seed(123)
index_train <- caret::createDataPartition(Y$event,
                                          p = .7, #70 %
                                          list = FALSE,
                                          times = 1)

X_train <- list()
X_test <- list()
for(omic in names(X)){
  X_train[[omic]] <- X[[omic]][index_train,,drop=F]
  X_test[[omic]] <- X[[omic]][-index_train,,drop=F]
}

Y_train <- Y[index_train,]
Y_test <- Y[-index_train,]
```

EVP per block:


```r
EPV <- getEPV.mb(X_train, Y_train)
for(b in names(X_train)){
  message(paste0("EPV = ", round(EPV[[b]], 4), ", for block ", b))
}
#> EPV = 0.0794, for block mirna
#> EPV = 0.1382, for block proteomic
```

## Cross Validation

In order to perform survival analysis with our high-dimensional data, we have implemented a series of methods that utilize techniques such as Cox Elastic Net, to select a lower number of features applying a penalty or partial least squares (PLS) methodology in order to reduce the dimensionality of the input data. 

To evaluate the performance of these methods, we have implemented cross-validation, which allows us to estimate the optimal parameters for future predictions based on prediction metrics such as: **AIC**, **C-INDEX**, **I.BRIER** and **AUC**. By default AUC metric (Area under the ROC curve) is used with the "cenROC" evaluator as it has provided the best results in our tests. However, multiple **AUC** evaluators could be used: "risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C" and "smoothROCtime_I". 

Furthermore, a mix of multiple metrics could be used to obtain the optimal model. The user has to establish different weight for each metric and all of them will be consider to compute the optimal model (the total weight must sum 1).

In addition, we have included options for normalizing data, filtering variables, and setting the minimum EPV, as well as specific parameters for each method, such as the alpha value for Cox Elastic Net and the number of components for PLS models. Overall, our cross-validation methodology allows us to effectively analyze high-dimensional survival data and optimize our model parameters.

## Preprocessing

As classical and PLS approach can not be perform in multiomic data, several multiomic methods have been designed. The methods are divided into two categories:

1. **Single-Block approach**: Each omic will be studied individually to reduce their dimensionality in order to explain the survival. Then, integration of all omics will be run to obtain the final model.
2. **Multi-Omic-Block approach**: All omics will be studied together by multiblock sPLS functions. After reduce the dimensionality, the latent variables will be used to obtain the final model.

But first, we establish the center and scale for each block:


```r
x.center = c(mirna = T, proteomic = T) #if vector, must be named
x.scale = c(mirna = F, proteomic = F) #if vector, must be named
```

## SB.sPLS-ICOX

In our pursuit of addressing the challenges posed by high-dimensional multi-omic data, we have developed SB.sPLS-ICOX, an innovative algorithm that combines the strengths of sPLS-ICOX and integrative analysis. SB.sPLS-ICOX employs a single-block approach, applying sPLS-ICOX individually to each omic, resulting in dimensionality reduction within each dataset. By constructing weights based on univariate cox models, we capture survival information during the reduction process. The reduced omic datasets are then integrated to create a comprehensive survival model using their PLS components. Cross-validation is utilized to determine the optimal number of components and the penalty for variable selection, ensuring robust model optimization.


```r
cv.sb.splsicox_res <- cv.sb.splsicox(X = X_train, Y = Y_train,
                                     max.ncomp = 2, penalty.list = c(0.5,0.9),
                                     n_run = 2, k_folds = 5, 
                                     x.center = x.center, x.scale = x.scale, 
                                     remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                                     remove_variance_at_fold_level = F,
                                     remove_non_significant_models = F, alpha = 0.05, 
                                     w_AIC = 0, w_c.index = 0, w_AUC = 1, w_BRIER = 0, times = NULL, max_time_points = 15,
                                     MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                                     pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
                                     MIN_EPV = 5, return_models = F, remove_non_significant = F, returnData = F,
                                     PARALLEL = F, verbose = F, seed = 123)

cv.sb.splsicox_res
```


```r
cv.sb.splsicox_res$plot_AUC
```

We will generate a SB.sPLS-ICOX model with optimal number of principal components and its penalty based on the results obtained from the cross validation.


```r
sb.splsicox_model <- sb.splsicox(X = X_train, Y = Y_train,
                                 n.comp = 1, #cv.sb.splsicox_res$opt.comp, 
                                 penalty = 0.9, #cv.sb.splsicox_res$opt.penalty,
                                 x.center = x.center, x.scale = x.scale,
                                 remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                                 remove_non_significant = F, 
                                 alpha = 0.05, MIN_EPV = 5, 
                                 returnData = T, verbose = F)
#> As we are working with a multiblock approach with 2 blocks, a maximum of 5 components can be used.

sb.splsicox_model
#> The method used is SB.sPLS-ICOX.
#> Survival model:
#>                        coef exp(coef)   se(coef)  robust se          z     Pr(>|z|)
#> comp_1_mirna      0.1077358 1.1137534 0.01865296 0.01661403  6.4846281 8.895092e-11
#> comp_1_proteomic -0.1416963 0.8678848 0.23679243 0.24770216 -0.5720431 5.672928e-01
```

In case some components get a P-Value greater than the cutoff for significant, we can drop them by the parameter "remove_non_significant".


```r
sb.splsicox_model <- sb.splsicox(X = X_train, Y = Y_train,
                                 n.comp = 1, #cv.sb.splsicox_res$opt.comp,
                                 penalty = 0.9, #cv.sb.splsicox_res$opt.penalty,
                                 x.center = x.center, x.scale = x.scale,
                                 remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                                 remove_non_significant = T,
                                 alpha = 0.05, MIN_EPV = 5, 
                                 returnData = T, verbose = F)
#> As we are working with a multiblock approach with 2 blocks, a maximum of 5 components can be used.

sb.splsicox_model
#> The method used is SB.sPLS-ICOX.
#> A total of 1 variables have been removed due to non-significance filter inside cox model.
#> Survival model:
#>                    coef exp(coef)   se(coef)  robust se        z     Pr(>|z|)
#> comp_1_mirna 0.09954254  1.104665 0.01257348 0.01573434 6.326453 2.508607e-10
```

In this case, we optimized each omic/block to use the same number of components. But there is another methodology that allow to select a different number of components per block call "isb.splsicox".


```r
cv.isb.splsicox_res <- cv.isb.splsicox(X = X_train, Y = Y_train,
                                     max.ncomp = 2, penalty.list = c(0.5,0.9),
                                     n_run = 2, k_folds = 5, 
                                     x.center = x.center, x.scale = x.scale, 
                                     remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                                     remove_variance_at_fold_level = F,
                                     remove_non_significant_models = F, alpha = 0.05, 
                                     w_AIC = 0, w_c.index = 0, w_AUC = 1, w_BRIER = 0, times = NULL, max_time_points = 15,
                                     MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                                     pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
                                     MIN_EPV = 5, return_models = F, remove_non_significant = F, returnData = F,
                                     PARALLEL = F, verbose = F, seed = 123)

cv.isb.splsicox_res
```


```r
cv.isb.splsicox_res$list_cv_spls_models$mirna$plot_AUC
cv.isb.splsicox_res$list_cv_spls_models$proteomic$plot_AUC
```


```r
isb.splsicox_model <- isb.splsicox(X = X_train, Y = Y_train, cv.isb = cv.isb.splsicox_res,
                                      x.center = x.center, x.scale = x.scale, 
                                      remove_near_zero_variance = TRUE, remove_zero_variance = TRUE, toKeep.zv = NULL,
                                   remove_non_significant = FALSE, alpha = 0.05,
                                   MIN_EPV = 5, returnData = TRUE, verbose = FALSE)

isb.splsicox_model
```

## SB.sPLS-DRCOX-Dynamic

We have developed also SB.sPLS-DRCOX, a comprehensive algorithm that employs the individual sPLS-DRCOX approach, executing the algorithm on each omic dataset independently to achieve dimensionality reduction. By integrating the resulting components, a unified survival model is constructed, capturing the collective information from all omics. Similar to SB.sPLS-ICOX, cross-validation is employed to identify the optimal number of components and penalty for variable selection, ensuring robust model optimization.


```r
cv.sb.splsdrcox_res <- cv.sb.splsdrcox(X = X_train, Y = Y_train,
                                       max.ncomp = 2, vector = NULL,
                                       n_run = 2, k_folds = 10, 
                                       x.center = x.center, x.scale = x.scale,
                                       #y.center = FALSE, y.scale = FALSE,
                                       remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                                       remove_variance_at_fold_level = F,
                                       remove_non_significant_models = F, alpha = 0.05, 
                                       w_AIC = 0, w_c.index = 0, w_AUC = 1, w_BRIER = 0, times = NULL, max_time_points = 15,
                                       MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                                       pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
                                       MIN_EPV = 5, return_models = F, remove_non_significant = F, returnData = F,
                                       PARALLEL = F, verbose = F, seed = 123)

cv.sb.splsdrcox_res
```

We will generate a SB.sPLS-DRCOX model with optimal number of principal components and its penalty based on the results obtained from the cross validation.


```r
sb.splsdrcox_model <- sb.splsdrcox(X = X_train, 
                                   Y = Y_train, 
                                   n.comp = 2, #cv.sb.splsdrcox_res$opt.comp, 
                                   vector = list("mirna" = 484, "proteomic" = 369), #cv.sb.splsdrcox_res$opt.nvar,
                                   x.center = x.center, x.scale = x.scale,
                                   remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL, 
                                   remove_non_significant = T, alpha = 0.05, MIN_EPV = 5,
                                   returnData = T, verbose = F)
#> As we are working with a multiblock approach with 2 blocks, a maximum of 5 components can be used.

sb.splsdrcox_model
#> The method used is SB.sPLS-DRCOX-Dynamic.
#> Survival model:
#>                          coef exp(coef)     se(coef)    robust se        z     Pr(>|z|)
#> comp_2_mirna     1.143874e-06  1.000001 5.309299e-07 4.324450e-07 2.645132 8.165901e-03
#> comp_1_proteomic 6.761452e-01  1.966283 9.738196e-02 1.015806e-01 6.656243 2.809162e-11
#> comp_2_proteomic 2.075853e-01  1.230703 5.953011e-02 5.226568e-02 3.971732 7.135210e-05
```

As before, we optimized each omic/block to use the same number of components. But there is another methodology that allow to select a different number of components and penalty per block call "cv.isb.splsdrcox".


```r
cv.isb.splsdrcox_res <- cv.isb.splsdrcox(X = X_train, Y = Y_train,
                                     max.ncomp = 2, vector = NULL,
                                     n_run = 2, k_folds = 5, 
                                     x.center = x.center, x.scale = x.scale, 
                                     remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                                     remove_variance_at_fold_level = F,
                                     remove_non_significant_models = F, alpha = 0.05, 
                                     w_AIC = 0, w_c.index = 0, w_AUC = 1, w_BRIER = 0, times = NULL, max_time_points = 15,
                                     MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                                     pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
                                     MIN_EPV = 5, return_models = F, remove_non_significant = F, returnData = F,
                                     PARALLEL = F, verbose = F, seed = 123)

cv.isb.splsdrcox_res
```


```r
cv.isb.splsdrcox_res$list_cv_spls_models$mirna$plot_AUC
cv.isb.splsdrcox_res$list_cv_spls_models$proteomic$plot_AUC
```


```r
isb.splsdrcox_model <- isb.splsdrcox(X = X_train, Y = Y_train, cv.isb = cv.isb.splsdrcox_res,
                                   x.center = x.center, x.scale = x.scale, 
                                   remove_near_zero_variance = TRUE, remove_zero_variance = TRUE, toKeep.zv = NULL,
                                   remove_non_significant = FALSE, alpha = 0.05,
                                   MIN_EPV = 5, returnData = TRUE, verbose = FALSE)

isb.splsdrcox_model
```

## MB.sPLS-DRCOX

To enhance the versatility and performance of our survival analysis methods, we have developed full multiblock survival models by combining the sPLS-DRCOX methodology with the multiblock sPLS functions from the mixOmics R package. This allows a full integration for omic by selecting the components to reach the same objective.

In the creation of these methods, we utilized an heuristic variable selection approach along with the MixOmics algorithms for hyperparameter selection. The penalty is determined based on a vector of the number of variables to test. Users have the flexibility to specify a specific value for selecting a fixed number of variables. Alternatively, by setting the value to NULL, the heuristic process automatically selects the optimal number of variables.

Our method empowers users to define the minimum and maximum number of variables to consider, as well as the number of cutpoints to test between these limits. Through an iterative process, the algorithm identifies the optimal number of variables and further explores the performance of existing cutpoints compared to the selected optimal value.

This integration of sPLS-DRCOX and multiblock sPLS provides researchers with a powerful tool for conducting comprehensive multivariate survival analysis.


```r
cv.mb.splsdrcox_res <- cv.mb.splsdrcox(X = X_train, Y = Y_train, 
                                       max.ncomp = 2, vector = NULL, #NULL - autodetection
                                       MIN_NVAR = 10, MAX_NVAR = NULL, n.cut_points = 10, EVAL_METHOD = "AUC",
                                       n_run = 2, k_folds = 4, 
                                       x.center = x.center, x.scale = x.scale, 
                                       remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                                       remove_variance_at_fold_level = F,
                                       remove_non_significant_models = F, alpha = 0.05, 
                                       w_AIC = 0, w_c.index = 0, w_AUC = 1, w_BRIER = 0, times = NULL, max_time_points = 15,
                                       MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                                       pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
                                       MIN_EPV = 5, return_models = F, remove_non_significant = F, returnData = F,
                                       PARALLEL = F, verbose = F, seed = 123)

cv.mb.splsdrcox_res
```

After getting the cross-validation, the full model could be obtain by passing the optimized number of components and the list of the number of variables per omic.


```r
mb.splsdrcox_model <- mb.splsdrcox(X = X_train, Y = Y_train, 
                                   n.comp = 2, #cv.mb.splsdrcox_res$opt.comp,
                                   vector = list("mirna" = 326, "proteomic" = 369), #cv.mb.splsdrcox_res$opt.nvar,
                                   x.center = x.center, x.scale = x.scale, 
                                   remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL,
                                   remove_non_significant = T, alpha = 0.05,
                                   MIN_AUC_INCREASE = 0.01,
                                   pred.method = "cenROC", max.iter = 200,
                                   times = NULL, max_time_points = 15,
                                   MIN_EPV = 5, returnData = T, verbose = F)
#> As we are working with a multiblock approach with 2 blocks, a maximum of 5 components can be used.
#> Design matrix has changed to include Y; each block will be
#>             linked to Y.

mb.splsdrcox_model
#> The method used is MB.sPLS-DRCOX-Dynamic.
#> A total of 3 variables have been removed due to non-significance filter inside cox model.
#> Survival model:
#>                      coef exp(coef)   se(coef)    robust se        z  Pr(>|z|)
#> comp_2_mirna 3.820833e-07         1 3.4797e-07 3.182173e-07 1.200699 0.2298678
```

In some cases, if any component is significant, the model will keep at least the one with the lower P-Value.

## MB.sPLS-DACOX

We also have extended our methodological repertoire with the development of MB.sPLS-DACOX. Building upon the foundation of sPLS-DACOX, this novel approach incorporates the powerful MultiBlock functions from the MixOmics R package, further enhancing its capabilities and performance.


```r
# run cv.splsdrcox
cv.mb.splsdacox_res <- cv.mb.splsdacox(X = X_train, Y = Y_train, 
                                       max.ncomp = 2, vector = NULL, #NULL - autodetection
                                       n_run = 2, k_folds = 4, 
                                       x.center = x.center, x.scale = x.scale, 
                                       remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                                       remove_variance_at_fold_level = F,
                                       remove_non_significant_models = F, alpha = 0.05, 
                                       w_AIC = 0, w_c.index = 0, w_AUC = 1, w_BRIER = 0, times = NULL, max_time_points = 15,
                                       MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                                       pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
                                       MIN_EPV = 5, return_models = F, remove_non_significant = F, returnData = F,
                                       PARALLEL = F, verbose = F, seed = 123)

cv.mb.splsdacox_res
```

After getting the cross-validation, the full model could be obtain by passing the optimized number of components and the list of the number of variables per omic.


```r
mb.splsdacox_model <- mb.splsdacox(X = X_train, Y = Y_train, 
                                   n.comp = 2, #cv.mb.splsdacox_res$opt.comp,
                                   vector = list("mirna" = 326, "proteomic" = 10), #cv.mb.splsdacox_res$opt.nvar,
                                   x.center = x.center, x.scale = x.scale, 
                                   remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL,
                                   remove_non_significant = T, alpha = 0.05,
                                   MIN_AUC_INCREASE = 0.01,
                                   pred.method = "cenROC", max.iter = 200,
                                   times = NULL, max_time_points = 15,
                                   MIN_EPV = 5, returnData = T, verbose = F)
#> As we are working with a multiblock approach with 2 blocks, a maximum of 5 components can be used.
#> Design matrix has changed to include Y; each block will be
#>             linked to Y.

mb.splsdacox_model
#> The method used is MB.sPLS-DACOX-Dynamic.
#> A total of 3 variables have been removed due to non-significance filter inside cox model.
#> Survival model:
#>                        coef exp(coef)   se(coef)  robust se        z  Pr(>|z|)
#> comp_2_proteomic 0.08820783  1.092215 0.07485843 0.07740266 1.139597 0.2544543
```

# Comparing multiple survival models

Next, we will analyze the results obtained from the multiple models to see which one obtains the best predictions based on our data. To do this, we will use the test set that has not been used for the training of any model.

## Comparing for multiple evaluators at the same time.

Initially, we will compare the area under the curve (AUC) for each of the methods according to the evaluator we want. The function is developed to simultaneously evaluate multiple evaluators. However, we will continue working with a single evaluator. In this case "cenROC". On the other hand, we must provide a list of the different models as well as the X and Y test set we want to evaluate.

When evaluating survival model results, we must indicate at which temporal points we want to perform the evaluation. As we already specified a NULL value for the "times" variable in the cross-validation, we are going to let the algorithm to compute them again.


```r
lst_models <- list("SB.sPLS-ICOX" = sb.splsicox_model,
                   #"iSB.sPLS-ICOX" = isb.splsicox_model,
                   "SB.sPLS-DRCOX-Dynamic" = sb.splsdrcox_model,
                   #"iSB.sPLS-DRCOX-Dynamic" = isb.splsdrcox_model,
                   #"SB.sPLS-DRCOX-Penalty" = sb.splsdrcox_penalty_model,
                   #"iSB.sPLS-DRCOX-Penalty" = isb.splsdrcox_penalty_model,
                   #"SB.sPLS-DACOX-Dynamic" = sb.splsdacox_model,
                   #"iSB.sPLS-DACOX-Dynamic" = isb.splsdacox_model,
                   "MB.sPLS-DRCOX" = mb.splsdrcox_model,
                   "MB.sPLS-DACOX" = mb.splsdacox_model)

eval_results <- eval_Coxmos_models(lst_models = lst_models,
                                  X_test = X_test, Y_test = Y_test, 
                                  pred.method = "cenROC",
                                  pred.attr = "mean",
                                  times = NULL, max_time_points = 15, 
                                  PARALLEL = F)
```

In case you prefer to test multiple AUC evaluators, a list of them could be proportionate by using the library "purrr".


```r
lst_evaluators <- c(cenROC = "cenROC", risksetROC = "risksetROC")

eval_results <- purrr::map(lst_evaluators, ~eval_Coxmos_models(lst_models = lst_models,
                                                              X_test = X_test, Y_test = Y_test,
                                                              pred.method = .,
                                                              pred.attr = "mean",
                                                              times = NULL, 
                                                              max_time_points = 15,
                                                              PARALLEL = F))
```

We can print the results obtained in the console where we can see, for each of the selected methods, the training time and the time it took to be evaluated, as well as their AIC, C-Index and AUC metrics and at which time points it was evaluated.


```r
eval_results
#> Evaluation performed for methods: SB.sPLS-ICOX, SB.sPLS-DRCOX-Dynamic, MB.sPLS-DRCOX, MB.sPLS-DACOX.
#> SB.sPLS-ICOX:
#> 	training.time: 0.1211
#> 	evaluating.time: 0.0101
#> 	AIC: 310.0196
#> 	c.index: 0.8146
#> 	time: 116, 369.786, 623.572, 877.358, 1131.143, 1384.929, 1638.715, 1892.5, 2146.286, 2400.072, 2653.858, 2907.643, 3161.429, 3415.215, 3669
#> 	AUC: 0.52672
#> 	brier_time: 116, 186, 322, 377, 385, 454, 477, 489, 506, 518, 522, 606, 614, 616, 678, 703, 723, 747, 759, 825, 847, 904, 1015, 1048, 1150, 1174, 1189, 1233, 1275, 1324, 1508, 1528, 1611, 1694, 1742, 1793, 1935, 1972, 2632, 2854, 2866, 3001, 3472, 3669
#> 	I.Brier: 0.18894
#> 
#> SB.sPLS-DRCOX-Dynamic:
#> 	training.time: 0.0436
#> 	evaluating.time: 0.0082
#> 	AIC: 318.4197
#> 	c.index: 0.8215
#> 	time: 116, 369.786, 623.572, 877.358, 1131.143, 1384.929, 1638.715, 1892.5, 2146.286, 2400.072, 2653.858, 2907.643, 3161.429, 3415.215, 3669
#> 	AUC: 0.60663
#> 	brier_time: 116, 186, 322, 377, 385, 454, 477, 489, 506, 518, 522, 606, 614, 616, 678, 703, 723, 747, 759, 825, 847, 904, 1015, 1048, 1150, 1174, 1189, 1233, 1275, 1324, 1508, 1528, 1611, 1694, 1742, 1793, 1935, 1972, 2632, 2854, 2866, 3001, 3472, 3669
#> 	I.Brier: 0.18619
#> 
#> MB.sPLS-DRCOX:
#> 	training.time: 0.0313
#> 	evaluating.time: 0.009
#> 	AIC: 368.5598
#> 	c.index: 0.5319
#> 	time: 116, 369.786, 623.572, 877.358, 1131.143, 1384.929, 1638.715, 1892.5, 2146.286, 2400.072, 2653.858, 2907.643, 3161.429, 3415.215, 3669
#> 	AUC: 0.47539
#> 	brier_time: 116, 186, 322, 377, 385, 454, 477, 489, 506, 518, 522, 606, 614, 616, 678, 703, 723, 747, 759, 825, 847, 904, 1015, 1048, 1150, 1174, 1189, 1233, 1275, 1324, 1508, 1528, 1611, 1694, 1742, 1793, 1935, 1972, 2632, 2854, 2866, 3001, 3472, 3669
#> 	I.Brier: 0.17287
#> 
#> MB.sPLS-DACOX:
#> 	training.time: 0.0602
#> 	evaluating.time: 0.0195
#> 	AIC: 368.8973
#> 	c.index: 0.5307
#> 	time: 116, 369.786, 623.572, 877.358, 1131.143, 1384.929, 1638.715, 1892.5, 2146.286, 2400.072, 2653.858, 2907.643, 3161.429, 3415.215, 3669
#> 	AUC: 0.57775
#> 	brier_time: 116, 186, 322, 377, 385, 454, 477, 489, 506, 518, 522, 606, 614, 616, 678, 703, 723, 747, 759, 825, 847, 904, 1015, 1048, 1150, 1174, 1189, 1233, 1275, 1324, 1508, 1528, 1611, 1694, 1742, 1793, 1935, 1972, 2632, 2854, 2866, 3001, 3472, 3669
#> 	I.Brier: 0.17014
#> 
```

## Plot comparison

However, we can also obtain graphical results where we can compare each method over time, as well as their average scores using the function "plot_evaluation" or "plot_evaluation.list" if multiple evaluators have been tested. The user could choose to plot the AUC for time prediction points or Brier Score. In case of use Brier Score, instead of uses the Integrative Brier Score for the boxplots, the mean or median is used (plot_evaluation parameter).


```r
lst_eval_results <- plot_evaluation(eval_results, evaluation = "AUC")
lst_eval_results_brier <- plot_evaluation(eval_results, evaluation = "Brier")
```

After performing the cross-validation, we obtain a list in R that contains two new lists. The first of these refers to the evaluation over time for each of the methods used, as well as a variant where the average value of each of them is shown. On the other hand, we can compare the mean results of each method using: T-test, Wilcoxon-test, anova or Kruskal-Wallis.


```r
lst_eval_results$lst_plots$lineplot.mean
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-171-1.png" width="100%" />

```r
lst_eval_results$lst_plot_comparisons$t.test
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-171-2.png" width="100%" />

```r

# lst_eval_results$cenROC$lst_plots$lineplot.mean
# lst_eval_results$cenROC$lst_plot_comparisons$t.test
```

## Computing time comparison (we use the cross validation models)

Another possible comparison is related to the computation times for cross-validation, model creation, and evaluation. In this case, cross-validations and methods are loaded.


```r
lst_models_time <- list(#cv.sb.splsicox_res,
                        sb.splsicox_model,
                        #isb.splsicox_model,
                        #cv.sb.splsdrcox_res,
                        sb.splsdrcox_model,
                        #isb.splsdrcox_model,
                        #cv.mb.splsdrcox_res,
                        mb.splsdrcox_model,
                        #cv.mb.splsdrcox_res,
                        mb.splsdacox_model,
                        eval_results)
```


```r
ggp_time <- plot_time.list(lst_models_time)
```


```r
ggp_time
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-174-1.png" width="100%" />

# Individual interpretations of the results

Following the cross validation, we have selected the sPLS-DACOX methodology as the most suitable model for our data. We will now study and interpret its results based on the study variables or latent variables used. In this case, we will examine some graphs of the model.

## Forest plots

A forest plot can be obtained as the first graph using the survminer R package. However, the function has been restructured to allow for simultaneous launch of an Coxmos class model or a list of Coxmos models using the `plot_forest()` or `plot_forest.list()` function.


```r
#lst_forest_plot <- plot_forest.list(lst_models)
lst_forest_plot <- plot_forest(lst_models$`SB.sPLS-DRCOX`)
```


```r
#lst_forest_plot$`SB.sPLS-DRCOX`
lst_forest_plot
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-176-1.png" width="100%" />

## PH Assumption

The following graph is related to one of the assumptions of the Cox models, called proportional hazard. 

In a Cox proportional hazards model, the proportional hazards assumption states that the hazard ratio (the risk of experiencing the event of interest) is constant over time for a given set of predictor variables. This means that the effect of the predictors on the hazard ratio does not change over time. This assumption is important because it allows for the interpretation of the model's coefficients as measures of the effect of the predictors on the hazard ratio. Violations of the proportional hazards assumption can occur if the effect of the predictors on the hazard ratio changes over time or if there is an interaction between the predictors and time. In these cases, the coefficients of the model may not accurately reflect the effect of the predictors on the hazard ratio and the results of the model may not be reliable.

In this way, to visualize and check if the assumption is violated, the function `plot_proportionalHazard.list()` or `plot_proportionalHazard()` can be called, depending on whether a list of models or a specific model is to be evaluated.


```r
#lst_ph_ggplot <- plot_proportionalHazard.list(lst_models)
lst_ph_ggplot <- plot_proportionalHazard(lst_models$`SB.sPLS-DRCOX`)
```

Variables or components with a significant P-Value indicate that the assumption is being violated.


```r
#lst_ph_ggplot$`SB.sPLS-DRCOX`
lst_ph_ggplot
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-178-1.png" width="100%" />

## Density plots

Another type of graph implemented for all models, whether they belong to the classical branch or to PLS-based models, is the visualization of observations by event according to the values predicted by the Cox models. 

The R package "coxph" allows for several types of predictions to be made on a Cox model that we use in our function, which are:

* Linear predictors "lp": are the expected values of the response variable (in this case, time until the event of interest) for each observation, based on the Cox model. These values can be calculated from the mean of the predictor variable values and the constant term of the model.

* Risk of experiencing an event "risk": is a measure of the probability that an event will occur for each observation, based on the Cox model. The risk value can be calculated from the predictor values and the constant term of the model.

* Number of events expected to be experienced over time with these specific individual characteristics "expected": are the expected number of events that would occur for each observation, based on the Cox model and a specified period of time.

* Terms: are the variables included in the Cox model.

* Survival probability "survival": is the probability that an individual will not experience the event of interest during a specified period of time, based on the Cox model. The survival probability can be calculated from the predictor values and the constant term of the model.

According to the predicted value, we can classify the observations along their possible values and see their distribution for each of the different models.


```r
#density.plots.lp <- plot_cox.event.list(lst_models, type = "lp")
density.plots.lp <- plot_cox.event(lst_models$`SB.sPLS-DRCOX`, type = "lp")
```


```r
density.plots.lp$plot.density
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-180-1.png" width="100%" />

```r
density.plots.lp$plot.histogram
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-180-2.png" width="100%" />

## Studying PLS model

For those models based on PLS components, the PLS could be studied in terms of loadings/scores. In order to get the plots, the function `plot_PLS_Coxmos()` has been developed where the user could specify to see "scores", "loadings" or a "biplot" for a couple of latent variables. By default, if no factor is given, samples are grouped by event.


```r
ggp_scores <- plot_PLS_Coxmos(model = lst_models$`SB.sPLS-DRCOX`,
                             comp = c(1,2), mode = "scores")
#> The model has only 1 component
```


```r
ggp_scores$plot_block
#> $mirna
#> $mirna$plot
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-182-1.png" width="100%" />

```
#> 
#> $mirna$outliers
#> NULL
#> 
#> 
#> $proteomic
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-182-2.png" width="100%" />


```r
ggp_loadings <- plot_PLS_Coxmos(model = lst_models$`SB.sPLS-DRCOX`, 
                               comp = c(1,2), mode = "loadings",
                               top = 10)
#> The model has only 1 component
```


```r
ggp_loadings$plot_block
#> $mirna
#> $mirna$plot
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-184-1.png" width="100%" />

```
#> 
#> $mirna$outliers
#> NULL
#> 
#> 
#> $proteomic
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-184-2.png" width="100%" />


```r
ggp_biplot <- plot_PLS_Coxmos(model = lst_models$`SB.sPLS-DRCOX`, 
                             comp = c(1,2), mode = "biplot",
                             top = 15,
                             only_top = T)
#> The model has only 1 component
```


```r
ggp_biplot$plot_block
#> $mirna
#> $mirna$plot
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-186-1.png" width="100%" />

```
#> 
#> $mirna$outliers
#> NULL
#> 
#> 
#> $proteomic
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-186-2.png" width="100%" />

# Understanding the results in terms of the original variables (PLS models)

When a PLS-COX model is computed, the final survival model is based on the PLS latent variables. Although those new variables can explain the survival, in order to understand the disease, new coefficients could be computed in terms of the original variables. 

Coxmos is able to transfer the component beta coefficient to original variables by decomposing the coefficients by using the weight of the variables in that latent variables.

## Understanding the relative contribution per variable/component

However, before studying the original variables, if a PLS model is computed. Coxmos also proportionates plots to see how many % of AUC is computed per each component in order to see which components or latent variables are more related to the observation survival.

The % of AUC explanation per component could be calculated for TRAIN or TEST data. Train data was used to compute the best model, so the sum of variables/components will generate a better LP (linear predictor) performance. However, when test data is used, could happen that a specific variable or a component could be a better predictor than the full model.


```r
variable_auc_results <- eval_Coxmos_model_per_variable(model = lst_models$`SB.sPLS-DRCOX`, 
                                                      X_test = lst_models$`SB.sPLS-DRCOX`$X_input, 
                                                      Y_test = lst_models$`SB.sPLS-DRCOX`$Y_input,
                                                      pred.method = "cenROC", pred.attr = "mean",
                                                      times = NULL, max_time_points = 15,
                                                      PARALLEL = FALSE)

variable_auc_plot_train <- plot_evaluation(variable_auc_results, evaluation = "AUC")
```


```r
variable_auc_plot_train$lst_plots$lineplot.mean
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-188-1.png" width="100%" />

The plot shows the AUC for the full model (called LP), and then, the AUC per each variable (for classical methods) or components (for PLS methods).

## Psudobeta

In order to improve the interpretability of a PLS model, a subset of the most influential variables can be selected. In this example, the top 20 variables are chosen. Additionally, non-significant PLS components are excluded by setting the "onlySig" parameter to "TRUE".


```r
# ggp.simulated_beta <- plot_pseudobpenalty.list(lst_models = lst_models, 
#                                            error.bar = T, onlySig = T, alpha = 0.05, 
#                                            zero.rm = T, auto.limits = T, top = 20,
#                                            show_percentage = T, size_percentage = 2, verbose = F)

ggp.simulated_beta <- plot_pseudobeta(model = lst_models$`SB.sPLS-DRCOX`, 
                                      error.bar = T, onlySig = T, alpha = 0.05, 
                                      zero.rm = T, auto.limits = T, top = 20,
                                      show_percentage = T, size_percentage = 2)
```

The iSB.sPLS-DRCOX model was computed using a total of 2 components. Although these components were classified as dangerous for the observations (based on coefficients greater than one), certain variables within the components may still have a protective effect, depending on their individual weights.

The following plot illustrates the pseudo-beta coefficient for the original variables in the iSB.sPLS-DRCOX model. As only the top 20 variables are shown, in some case a percentage of the total linear predictor total value could be shown. To view the complete model, all variables would need to be included in the plot by assigning the value NULL to the "top" parameter.


```r
ggp.simulated_beta$plot
#> $mirna
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-190-1.png" width="100%" />

```
#> 
#> $proteomic
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-190-2.png" width="100%" />

## Kaplan-Meier

For a more intuitive understanding of the model, the user can also employ the `getAutoKM.list()` or `getAutoKM()` functions to generate Kaplan-Meier curves. These functions allow the user to view the KM curve for the entire model, a specific component of a PLS model, or for individual variables.

### Full model

To run the full Kaplan-Meier model, the "type" parameter must be set to "LP" (linear predictors). This means that the linear predictor value for each patient will be used to divide the groups (in this case, the score value multiplied by the Cox coefficient). The "surv_cutpoint" function from the R package "survminer" is used to determine the optimal cut-point. Other parameters are not used in this specific method.


```r
# LST_KM_RES_LP <- getAutoKM.list(type = "LP",
#                                 lst_models = lst_models,
#                                 comp = 1:4,
#                                 top = 10,
#                                 ori_data = T,
#                                 BREAKTIME = NULL,
#                                 only_sig = T, alpha = 0.05)

LST_KM_RES_LP <- getAutoKM(type = "LP",
                           model = lst_models$`SB.sPLS-DRCOX`,
                           comp = 1:4,
                           top = 10,
                           ori_data = T,
                           BREAKTIME = NULL,
                           only_sig = T, alpha = 0.05)
```

As a result, the Kaplan-Meier curve could be plotted.  


```r
LST_KM_RES_LP$LST_PLOTS$LP
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-192-1.png" width="100%" />

After generating a Kaplan-Meier curve for the model, the cutoff value, which is used to divide the observations into two groups, can be used to evaluate how the test data is classified. The vector of cutoffs for multiple models, when `getAutoKM.list()` is applied, can be retrieved by using the function `getCutoffAutoKM.list()` and passing the output of `getAutoKM.list()` as a parameter.

Once the vector is obtained, the function `getTestKM.list()` or `getTestKM()` can be run with the list of models, the X test data, Y test data, the list of cutoffs or a single value, and the desired number of breaks for the new Kaplan-Meier plot.

A log-rank test will be displayed to determine if the chosen cutoff is an effective way to split the data into groups with higher and lower risk.


```r
# lst_cutoff <- getCutoffAutoKM.list(LST_KM_RES_LP)
# LST_KM_TEST_LP <- getTestKM.list(lst_models = lst_models, 
#                                  X_test = X_test, Y_test = Y_test, 
#                                  type = "LP",
#                                  BREAKTIME = NULL, n.breaks = 20,
#                                  lst_cutoff = lst_cutoff)

lst_cutoff <- getCutoffAutoKM(LST_KM_RES_LP)
LST_KM_TEST_LP <- getTestKM(model = lst_models$`SB.sPLS-DRCOX`, 
                            X_test = X_test, Y_test = Y_test, 
                            type = "LP",
                            BREAKTIME = NULL, n.breaks = 20,
                            cutoff = lst_cutoff)
```


```r
LST_KM_TEST_LP
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-194-1.png" width="100%" />

### Components

To generate a Kaplan-Meier curve for a specific component, the "type" parameter must be set to "COMP" (component). This means that the linear predictor is computed using only one component at a time to split the groups. In this case, the "comp" parameter can be used to specify which component should be computed (if multiple components, each one will be computed separately).


```r
# LST_KM_RES_COMP <- getAutoKM.list(type = "COMP",
#                                   lst_models = lst_models,
#                                   comp = 1:4,
#                                   top = 10,
#                                   ori_data = T,
#                                   BREAKTIME = NULL,
#                                   only_sig = T, alpha = 0.05)

LST_KM_RES_COMP <- getAutoKM(type = "COMP",
                             model = lst_models$`SB.sPLS-DRCOX`,
                             comp = 1:4,
                             top = 10,
                             ori_data = T,
                             BREAKTIME = NULL,
                             only_sig = T, alpha = 0.05)
```


```r
LST_KM_RES_COMP$LST_PLOTS$mirna$comp_2
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-196-1.png" width="100%" />

```r
LST_KM_RES_COMP$LST_PLOTS$proteomic$comp_1
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-196-2.png" width="100%" />

```r
LST_KM_RES_COMP$LST_PLOTS$proteomic$comp_2
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-196-3.png" width="100%" />


```r
# lst_cutoff <- getCutoffAutoKM.list(LST_KM_RES_COMP)
# LST_KM_TEST_COMP <- getTestKM.list(lst_models = lst_models, 
#                                    X_test = X_test, Y_test = Y_test, 
#                                    type = "COMP",
#                                    BREAKTIME = NULL, n.breaks = 20,
#                                    lst_cutoff = lst_cutoff)

lst_cutoff <- getCutoffAutoKM(LST_KM_RES_COMP)
LST_KM_TEST_COMP <- getTestKM(model = lst_models$`SB.sPLS-DRCOX`, 
                              X_test = X_test, Y_test = Y_test, 
                              type = "COMP",
                              BREAKTIME = NULL, n.breaks = 20,
                              cutoff = lst_cutoff)
```


```r
LST_KM_TEST_COMP$comp_2_mirna
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-198-1.png" width="100%" />

```r
LST_KM_TEST_COMP$comp_1_proteomic
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-198-2.png" width="100%" />

```r
LST_KM_TEST_COMP$comp_2_proteomic
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-198-3.png" width="100%" />

### Original variables

To generate a Kaplan-Meier curve for original variables, the "type" parameter must be set to "VAR" (variable). In this case, the "ori_data" parameter can be used to determine whether the original or normalized values of the variables should be used.

Additionally, the "top" parameter can be used to plot a specific number of variables, sorted by P-Value (Log-Rank test). The best cutpoint is determined by the "surv_cutpoint" function. Both qualitative and quantitative variables are supported.


```r
# LST_KM_RES_VAR <- getAutoKM.list(type = "VAR",
#                                  lst_models = lst_models,
#                                  comp = 1:4,
#                                  top = 10,
#                                  ori_data = T,
#                                  BREAKTIME = NULL,
#                                  only_sig = T, alpha = 0.05)

LST_KM_RES_VAR <- getAutoKM(type = "VAR",
                            model = lst_models$`SB.sPLS-DRCOX`,
                            comp = 1:4,
                            top = 10,
                            ori_data = T,
                            BREAKTIME = NULL,
                            only_sig = T, alpha = 0.05)
```


```r
LST_KM_RES_VAR$LST_PLOTS$mirna$hsa.minus.miR.minus.21.minus.5p
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-200-1.png" width="100%" />

```r
LST_KM_RES_VAR$LST_PLOTS$proteomic$var_840
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-200-2.png" width="100%" />

```r
LST_KM_RES_VAR$LST_PLOTS$proteomic$var_7535
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-200-3.png" width="100%" />


```r
# lst_cutoff <- getCutoffAutoKM.list(LST_KM_RES_VAR)
# LST_KM_TEST_VAR <- getTestKM.list(lst_models = lst_models, 
#                                   X_test = X_test, Y_test = Y_test, 
#                                   type = "VAR", ori_data = T,
#                                   BREAKTIME = NULL, n.breaks = 20,
#                                   lst_cutoff = lst_cutoff)

lst_cutoff <- getCutoffAutoKM(LST_KM_RES_VAR)
LST_KM_TEST_VAR <- getTestKM(model = lst_models$`SB.sPLS-DRCOX`, 
                             X_test = X_test, Y_test = Y_test, 
                             type = "VAR", ori_data = T,
                             BREAKTIME = NULL, n.breaks = 20,
                             cutoff = lst_cutoff)
```


```r
LST_KM_TEST_VAR$mirna$hsa.minus.miR.minus.21.minus.5p
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-202-1.png" width="100%" />

```r
LST_KM_TEST_VAR$proteomic$var_840
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-202-2.png" width="100%" />

```r
LST_KM_TEST_VAR$proteomic$var_7535
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-202-3.png" width="100%" />

# New patients

In addition, Coxmos can also manage new patients to perform predictions.

To demonstrate, an observation from a test dataset will be used.


```r
new_pat <- list()
for(b in names(X_test)){
  new_pat[[b]] <- X_test[[b]][1,,drop=F]
}
```

As shown, this is a censored patient who has the last observation at time 825.


```r
knitr::kable(Y_test[rownames(new_pat$mirna),])
```



|                 | time|event |
|:----------------|----:|:-----|
|TCGA-A2-A0SV-01A |  825|TRUE  |



The function `plot_pseudobeta_newObservation.list()` or `plot_pseudobeta_newObservation()` allows the user to compare the characteristics of a new patient with the pseudo-beta values obtained from a specific model. The goal is to understand in a general way which variables are associated with an increased risk of an event or a decreased risk, in comparison to the variables predicted by the model.


```r
# ggp.simulated_beta_newPat <- plot_pseudobeta_newObservation.list(lst_models = lst_models, 
#                                                              new_observation = new_pat,
#                                                              error.bar = T, onlySig = T, alpha = 0.05,
#                                                              zero.rm = T, auto.limits = T, show.betas = T, top = 20)

ggp.simulated_beta_newPat <- plot_pseudobeta_newObservation(model = lst_models$`SB.sPLS-DRCOX`,
                                                        new_observation = new_pat,
                                                        error.bar = T, onlySig = T, alpha = 0.05,
                                                        zero.rm = T, auto.limits = T, show.betas = T, top = 20)
```

On the left, a linear predictor value is shown for the observation and each variable. On the right, the pseudo-beta coefficients for the model's original variables are illustrated. This allows the user to compare the direction of the linear predictors per variable. A change in direction means that the variable's value is below the mean for hazard variables or above the mean for protective ones. Having an opposite direction for hazard variables and maintaining the direction for protective variables makes the observation safer over time.


```r
ggp.simulated_beta_newPat$plot$mirna
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-206-1.png" width="100%" />

```r
ggp.simulated_beta_newPat$plot$proteomic
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-206-2.png" width="100%" />

## Add patient to density plot

Additionally, the new observation can be added to the density and histogram plots that the model has computed. While this function can be useful in visualizing the results, it is important to note that no definitive conclusions can be drawn from these types of plots alone. They serve as an additional means of viewing the predictive values for the new patient ("lp", "risk", "expected", "terms", "survival") in relation to the training dataset and its histogram or density plot, along with the events or censored patients.


```r
pat_density <- plot_observation.eventDensity(observation = new_pat, 
                                             model = lst_models$`SB.sPLS-DRCOX`, 
                                             time = NULL, 
                                             type = "lp")
```


```r
pat_density
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-208-1.png" width="100%" />


```r
pat_histogram <- plot_observation.eventHistogram(observation = new_pat, 
                                                 model = lst_models$`SB.sPLS-DRCOX`, 
                                                 time = NULL, 
                                                 type = "lp")
```


```r
pat_histogram
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-210-1.png" width="100%" />

## COX compare patients

Furthermore, similar plots can be generated for multiple patients at the same time. For example, by selecting 5 patients, individual linear predictors for each variable and the final linear predictor using all the variables used in the survival model can be plotted.

As an idea, this function could also be used to study the same patient in different time points after receiving a treatment and study its own evolution related to the original variables.


```r
sub_X_test <- list()
for(b in names(X_test)){
  sub_X_test[[b]] <- X_test[[b]][1:5,]
}
```


```r
knitr::kable(Y_test[rownames(sub_X_test$proteomic),])
```



|                 | time|event |
|:----------------|----:|:-----|
|TCGA-A2-A0SV-01A |  825|TRUE  |
|TCGA-A2-A0YT-01A |  723|TRUE  |
|TCGA-BH-A1FB-01A | 3669|TRUE  |
|TCGA-D8-A1XC-01A |  377|TRUE  |
|TCGA-BH-A1EX-01A | 1508|TRUE  |



To generate the figure, the function `plot_LP.multipleObservations.list()` or `plot_LP.multipleObservations()` can be used.


```r
# lst_cox.comparison <- plot_LP.multipleObservations.list(lst_models = lst_models,
#                                                     new_observations = sub_X_test,
#                                                     error.bar = T, zero.rm = T, onlySig = T,
#                                                     alpha = 0.05, top = 5)

lst_cox.comparison <- plot_LP.multipleObservations(model = lst_models$`SB.sPLS-DRCOX`,
                                                   new_observations = sub_X_test,
                                                   error.bar = T, zero.rm = T, onlySig = T,
                                                   alpha = 0.05, top = 5)
```


```r
lst_cox.comparison$plot
#> $mirna
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-214-1.png" width="100%" />

```
#> 
#> $proteomic
```

<img src="Coxmos-MO-pipeline_files/figure-html/unnamed-chunk-214-2.png" width="100%" />
