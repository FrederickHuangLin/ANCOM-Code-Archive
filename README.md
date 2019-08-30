# User Manual for ANCOM v2.0

The current code implements ANCOM in cross-sectional and longitudinal datasets while allowing the use of covariates. The following libraries need to be included for the R code to run:

```r
library(exactRankTests)
library(nlme)
library(dplyr)
```

## Input requirements

1. OTU data or taxa data: This is should be a matrix or data frame with each sample in rows and OTUs in columns. OTU data should contains a sample identifier with column name “```Sample.ID```”.
2. Metadata: This is the file with all variables and covariates of interest. It should be a matrix or data frame containing a sample identifier named “```Sample.ID```” and each following column being the variables.

## Arguments of the function (ANCOM)

* ```otu_data```: Data frame representing OTU data or taxa data (must meet specifications as mentioned earlier).
* ```meta_data```: Data frame of variables (must meet specifications as mentioned earlier).
* ```main_var```: Character. The name of the main variable of interest. 
* ```zero_cut```: Numeric fraction. OTUs with proportion of zeroes greater than prev.cut are not included in the analysis. Default is 0.90.
* ```p_adjust_method```: Character. Specifying the method to adjust p-values for multiple comparisons. Default is “BH” (Benjamini-Hochberg procedure).
* ```alpha```: Level of significance. Default is 0.05.
* ```adj_formula```: Character string representing the formula for adjustment (see example).
* ```rand_formula```: Character string representing the formula for random effects in lme (see example).



