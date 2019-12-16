# User Manual for [ANCOM](https://www.tandfonline.com/doi/full/10.3402/mehd.v26.27663) v2.0

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3472100.svg)](https://doi.org/10.5281/zenodo.3472100)

The current code implements ANCOM in cross-sectional and longitudinal datasets while allowing the use of covariates. The following libraries need to be included for the R code to run:

```r
library(exactRankTests)
library(nlme)
library(dplyr)
```

## Instructions for use

### Data preprocess

#### Usage

* ```feature_table_pre_process(feature_table, meta_data, sample_var, group_var = NULL, out_cut = 0.05, zero_cut = 0.90, lib_cut = 1000, neg_lb)```

#### Arguments

*	```feature_table```: Data frame or matrix representing observed OTU table with OTUs (or taxa) in rows and samples in columns.
*	```meta_data```: Data frame or matrix of all variables and covariates of interest.
*	```sample_var```: Character. The name of column storing sample IDs.
*	```group_var```: Character. The name of the group indicator. ```group_var``` is required for detecting structural zeros and outliers.
*	```out_cut``` Numerical fraction between 0 and 1. For each taxon, observations with proportion of mixture distribution less than ```out_cut``` will be detected as outlier zeros; while observations with proportion of mixture distribution greater than ```1 - out_cut``` will be detected as outlier values.
*	```zero_cut```: Numerical fraction between 0 and 1. Taxa with proportion of zeroes greater than ```zero_cut``` are not included in the analysis.
* ```lib_cut```: Numeric. Samples with library size less than ```lib_cut``` are not included in the analysis.
*	```neg_lb```: Logical. TRUE indicates a taxon would be classified as a structural zero in the corresponding experimental group using its asymptotic lower bound.

#### Value

* ```feature_table```: A data frame of pre-processed OTU table.
*	```meta_data```: A data frame of pre-processed metadata.
*	```structure_zeros```: A matrix consists of 0 and 1s with 1 indicating the taxon is identified as a structural zero in the corresponding group.

### ANCOM main function

#### Usage

* ```ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, alpha, adj_formula, rand_formula)```

####Arguments

* ```feature_table```: Data frame representing OTU/taxa data with OTUs (or taxa) in rows and samples in columns. Can be the output value from ```feature_table_pre_process```.
* ```meta_data```: Data frame of variables. Can be the output value from ```feature_table_pre_process```.
* ```struc_zero```: A matrix consists of 0 and 1s with 1 indicating the taxon is identified as a structural zero in the corresponding group. Can be the output value from ```feature_table_pre_process```.
* ```main_var```: Character. The name of the main variable of interest. 
* ```p_adjust_method```: Character. Specifying the method to adjust p-values for multiple comparisons. Default is “BH” (Benjamini-Hochberg procedure).
* ```alpha```: Level of significance. Default is 0.05.
* ```adj_formula```: Character string representing the formula for adjustment (see example).
* ```rand_formula```: Character string representing the formula for random effects in ```lme``` (see example).

#### A flowchart of the tests within ANCOM
![Flow Chart](/images/flowchart.png)

####Value

* ```res```: A data frame with the ```W``` statistic for each taxa and subsequent columns which are logical indicators of whether an OTU or taxon is differentially abundant under a series of cutoffs (0.9, 0.8, 0.7 and 0.6).

## Examples

### Standard analysis

_Detection of differentially abundant OTU between subjects_ <br/>
_Example dataset: moving-pics_

```r
# Step 1: Data preprocessing

feature_table = otu_data; sample_var = "Sample.ID"; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM

main_var = "Subject"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = NULL
t_start = Sys.time()
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start # around 30s

write_csv(res, "outputs/res_moving_pics.csv")

# Step 3 (Optional): Volcano Plot

# Calculate clr
clr_table = apply(feature_table, 2, clr)
# Calculate clr mean difference
eff_size = apply(clr_table, 1, function(y) 
  lm(y ~ x, data = data.frame(y = y, x = meta_data %>% pull(main_var)))$coef[-1])
# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Data frame for the figure
dat_fig = data.frame(taxa_id = res$taxa_id, x = eff_size, y = res$W) %>% 
  mutate(zero_ind = factor(ifelse(is.infinite(y), "Yes", "No"), levels = c("Yes", "No")))
# Replace Inf by (n_taxa - 1) for structural zeros
dat_fig$y = replace(dat_fig$y, is.infinite(dat_fig$y), n_taxa - 1)
# Annotation data
dat_ann = data.frame(x = min(dat_fig$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig = ggplot(data = dat_fig) + aes(x = x, y = y) + 
  geom_point(aes(color = zero_ind)) + 
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE) + 
  labs(x = "CLR mean difference", y = "W statistic", title = "Volcano Plot") +
  scale_color_discrete(name = "Structural\nzero", drop = FALSE) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
fig  
```

### Adjusted for covariates

_Detection of differentially abundant OTU between subjects adjusted for antibiotic usage_ <br/>
_Example dataset: moving-pics_

```r
main_var = "Subject"; p_adj_method = "BH"; alpha = 0.05
adj_formula = ”ReportedAntibioticUsage”; rand_formula = NULL
t_start = Sys.time()
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start
```

### Repeated measure/longitudinal analysis

#### Random intercept model

_Detection of differentially abundant genera between delivery methods accounting for random subject effect_ <br/>
_Each subject has his/her own intercept_ <br/>
_Example dataset: ecam_

```r
# Step 1: Data preprocessing

feature_table = otu_data; sample_var = "Sample.ID"; group_var = "delivery"
out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = TRUE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM

main_var = "delivery"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = "~ 1 | studyid"
t_start = Sys.time()
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start # around 30s

write_csv(res, "outputs/res_ecam.csv")

# Step 3 (Optional): Volcano Plot

# Calculate clr
clr_table = apply(feature_table, 2, clr)
# Calculate clr mean difference
eff_size = apply(clr_table, 1, function(y) 
  lm(y ~ x, data = data.frame(y = y, x = meta_data %>% pull(main_var)))$coef[-1])
# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Data frame for the figure
dat_fig = data.frame(taxa_id = res$taxa_id, x = eff_size, y = res$W) %>% 
  mutate(zero_ind = factor(ifelse(is.infinite(y), "Yes", "No"), levels = c("Yes", "No")))
# Replace Inf by (n_taxa - 1) for structural zeros
dat_fig$y = replace(dat_fig$y, is.infinite(dat_fig$y), n_taxa - 1)
# Annotation data
dat_ann = data.frame(x = min(dat_fig$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig = ggplot(data = dat_fig) + aes(x = x, y = y) + 
  geom_point(aes(color = zero_ind)) + 
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE) + 
  labs(x = "CLR mean difference", y = "W statistic", title = "Volcano Plot") +
  scale_color_discrete(name = "Structural\nzero", drop = FALSE) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
fig  
```

#### Random intercept model adjusted for other covariates

_Detection of differentially abundant genera between delivery methods accounting for fixed time effect and random subject effect_ <br/>
_Each subject has his/her own intercept_ <br/>
_Example dataset: ecam_

```r
main_var = "delivery"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "month"; rand_formula = "~ 1 | studyid"
t_start = Sys.time()
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start
```

#### Random coefficients/slope model

_Detection of differentially abundant genera between delivery methods accounting for random time effect and random subject effect_ <br/> 
_Each subject has his/her own intercept and slope_ <br/> 
_Example dataset: ecam_

```r
main_var = "delivery"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "month"; rand_formula = "~ month | studyid"
t_start = Sys.time()
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start
```






