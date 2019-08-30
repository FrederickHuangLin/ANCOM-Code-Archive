# User Manual for ANCOM v2.0

The current code implements ANCOM in cross-sectional and longitudinal datasets while allowing the use of covariates. The following libraries need to be included for the R code to run:

```r
library(exactRankTests)
library(nlme)
library(dplyr)
```

## Input requirements

1. OTU data or taxa data: This is should be a matrix or data frame with each sample in rows and OTUs in columns. OTU data should contains sample identifier with column name “```Sample.ID```”.
2. Metadata: This is the datafile with all variables and covariates of interest. It should be a matrix or data frame containing the sample identifier with column name “```Sample.ID```” and each following column being the variables.


