# User Manual for [ANCOM](https://www.tandfonline.com/doi/full/10.3402/mehd.v26.27663) v2.0

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3472100.svg)](https://doi.org/10.5281/zenodo.3472100)

The current code implements ANCOM in cross-sectional and longitudinal datasets while allowing the use of covariates. The following libraries need to be included for the R code to run:

```r
library(exactRankTests)
library(nlme)
library(dplyr)
```

## Input requirements

* OTU data or taxa data: This is should be a matrix or data frame with each sample in rows and OTUs in columns. OTU data should contains a sample identifier with column name “```Sample.ID```”.
* Metadata: This is the file with all variables and covariates of interest. It should be a matrix or data frame containing a sample identifier named “```Sample.ID```” and each following column being the variables.

## Arguments of the function (ANCOM)

* ```otu_data```: Data frame representing OTU data or taxa data (must meet specifications as mentioned earlier).
* ```meta_data```: Data frame of variables (must meet specifications as mentioned earlier).
* ```main_var```: Character. The name of the main variable of interest. 
* ```zero_cut```: Numeric fraction. OTUs with proportion of zeroes greater than ```zero_cut``` are not included in the analysis. Default is 0.90.
* ```p_adjust_method```: Character. Specifying the method to adjust p-values for multiple comparisons. Default is “BH” (Benjamini-Hochberg procedure).
* ```alpha```: Level of significance. Default is 0.05.
* ```adj_formula```: Character string representing the formula for adjustment (see example).
* ```rand_formula```: Character string representing the formula for random effects in ```lme``` (see example).

## A flowchart of the tests within ANCOM

![Flow Chart](/images/flowchart.png)

## Function outputs

* ```res```: A data frame with the ```W``` statistic for each taxa and subsequent columns which are logical indicators of whether an OTU or taxon is differentially abundant under a series of cutoffs (0.9, 0.8, 0.7 and 0.6).

## Examples

### Standard analysis

_Detection of differentially abundant OTU between subjects_ <br/>
_Example dataset: moving-pics_

```r
res = ANCOM(otu_data = otu_data, 
            meta_data = meta_data, 
            main_var = “Subject”,  
            zero_cut = 0.90, 
            p_adjust_method = “BH”, 
            alpha = 0.05, 
            adj_formula = NULL, 
            rand_formula = NULL)
```

### Adjusted for covariates

_Detection of differentially abundant OTU between subjects adjusted for antibiotic usage_ <br/>
_Example dataset: moving-pics_

```r
res = ANCOM(otu_data = otu_data, 
            meta_data = meta_data, 
            main_var = “Subject”,  
            zero_cut = 0.90, 
            p_adjust_method = “BH”, 
            alpha = 0.05, 
            adj_formula = ”ReportedAntibioticUsage”,
            rand_formula = NULL)
```

### Repeated measure/longitudinal analysis

#### Random intercept model

_Detection of differentially abundant genera between delivery methods accounting for random subject effect_ <br/>
_Each subject has his/her own intercept_ <br/>
_Example dataset: ecam_

```r
res = ANCOM(otu_data = otu_data, 
            meta_data = meta_data, 
            main_var = "delivery",  
            zero_cut = 0.90, 
            p_adjust_method = “BH”, 
            alpha = 0.05, 
            adj_formula = NULL,
            rand_formula = "~ 1 | studyid")
```

#### Random intercept model adjusted for other covariates

_Detection of differentially abundant genera between delivery methods accounting for fixed time effect and random subject effect_ <br/>
_Each subject has his/her own intercept_ <br/>
_Example dataset: ecam_

```r
res = ANCOM(otu_data = otu_data, 
            meta_data = meta_data, 
            main_var = "delivery",  
            zero_cut = 0.90, 
            p_adjust_method = “BH”, 
            alpha = 0.05, 
            adj_formula = ”month”,
            rand_formula = "~ 1 | studyid")
```

#### Random coefficients/slope model

_Detection of differentially abundant genera between delivery methods accounting for random time effect and random subject effect_ <br/> 
_Each subject has his/her own intercept and slope_ <br/> 
_Example dataset: ecam_

```r
res = ANCOM(otu_data = otu_data, 
            meta_data = meta_data, 
            main_var = "delivery",  
            zero_cut = 0.90, 
            p_adjust_method = “BH”, 
            alpha = 0.05, 
            adj_formula = ”month”,
            rand_formula = "~ month | studyid")
```






